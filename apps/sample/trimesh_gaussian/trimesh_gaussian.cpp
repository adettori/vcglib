#include <iostream>
#include <cmath>
#include <vcg/math/base.h>
#include <vcg/math/quaternion.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>
#include <wrap/io_trimesh/io_ply.h>
#include "./import_ply_GS.h"
#include "./export_ply_GS.h"
#include "./import_splat.h"
#include "export_splat.h"
#include "import_ksplat.h"

using namespace std;
using namespace vcg;

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>   ::AsVertexType,
                                           Use<MyEdge>     ::AsEdgeType,
                                           Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Color4b, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face<   MyUsedTypes, face::FFAdj,  face::Color4b, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<   MyUsedTypes> {};

class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};



int main(int argc, char *argv[])
{
    if(argc < 7) {
        cout << "Insufficient arguments" << endl;
        cout << "Expected args: input.ply minXBox minYBox minZBox maxXBox maxYBox maxZBox [theta] [phi]" << endl;
        return -1;
    }

    MyMesh mEllips, mSampledEllips;

    tri::io::PlyInfo pi;
    Matrix44f transf = Matrix44<float>();
    Quaternion<float> rotQuat;
    Point3<float> ellipseScale;
    Color4b ellipseColor;
    float theta = 0;
    float phi = 0;
    if(argc>9)
    {
        theta=stof(argv[8]);
        phi=stof(argv[9]);
    }

    MyMesh gauss;
    const int DegreeSH = 3;
    int ret = tri::io::ImporterKSPLAT<MyMesh, DegreeSH>::Open(gauss, argv[1]);
    if(ret != 0) {
        cout << "Error encountered while importing Gaussian splats: " << ret << endl;
    }

    MyMesh::PerVertexAttributeHandle<GaussianSplat<float,DegreeSH>> handleGauss =
        vcg::tri::Allocator<MyMesh>:: template GetPerVertexAttribute<GaussianSplat<float,DegreeSH>> (gauss, "gs");

    if(!tri::Allocator<MyMesh>::IsValidHandle(gauss, handleGauss))
    {
        return -1;
    }

    cout << "Number of Gaussian Splats: " << gauss.VN() << endl;
    Box3f box;
    if(argc>7)
    {
        box.min = Point3<float>(stof(argv[2]), stof(argv[3]), stof(argv[4]));
        box.max = Point3<float>(stof(argv[5]), stof(argv[6]), stof(argv[7]));
    }
    else // if not specified use the bbox of the mesh  to init the box
    {
        tri::UpdateBounding<MyMesh>::Box(gauss);
        Point3f c = gauss.bbox.Center();
        box.min = c - (c-gauss.bbox.min)*0.5;
        box.max = c + (gauss.bbox.max-c)*0.5;
    }

    // Loop through vertices, where each one represents a gaussian
    for(MyMesh::VertexIterator gi=gauss.vert.begin();gi!=gauss.vert.end();++gi) {

        if(tri::Index(gauss,*gi) % 10000 == 0) {
            printf("GS id: %9zu on %9zu (Deleted %8zu) \r",tri::Index(gauss,*gi), gauss.vert.size(), gauss.vert.size() - gauss.vn );
            fflush(stdout);
        }

        if(!box.IsIn(gi->P()) || GaussianSplat<float,DegreeSH>::getColor(handleGauss[gi])[3]<50) {
                tri::Allocator<MyMesh>::DeleteVertex(gauss, *gi);
        }
        else
        {
            // Splat is inside the box, continue the processing
            // Create solid
            tri::Octahedron(mEllips);

            tri::UpdatePosition<MyMesh>::Scale(mEllips, GaussianSplat<float,DegreeSH>::getScale(handleGauss[gi]));
            transf.SetIdentity();
            QuaternionToMatrix<float,Matrix44<float>>(GaussianSplat<float,DegreeSH>::getRotation(handleGauss[gi]), transf);
            tri::UpdatePosition<MyMesh>::Matrix(mEllips, transf);
            tri::UpdatePosition<MyMesh>::Translate(mEllips, gi->P());
            // Color
            tri::UpdateColor<MyMesh>::PerFaceConstant(mEllips, GaussianSplat<float,DegreeSH>::getColor(handleGauss[gi], theta, phi));
            tri::Append<MyMesh, MyMesh>::Mesh(mSampledEllips, mEllips);
            mEllips.Clear();
        }
    }

    tri::Allocator<MyMesh>::CompactVertexVector(gauss);
    printf("\n");

    printf("Saving surviving vertices (%i) \n", gauss.vn);
    tri::io::ExporterPLYGS<MyMesh, DegreeSH>::Save(gauss, "filteredGS.ply", true);

    printf("Saving Colored Ellipsoids\n");
    tri::io::ExporterPLY<MyMesh>::Save(mSampledEllips, "ellipseSamplesSH.ply", tri::io::Mask::IOM_FACECOLOR);

    return 0;
}