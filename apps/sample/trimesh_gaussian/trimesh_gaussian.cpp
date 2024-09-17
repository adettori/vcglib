#include <iostream>
#include <cmath>
#include <vcg/math/base.h>
#include <vcg/math/quaternion.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>
#include <wrap/io_trimesh/io_ply.h>
#include "./import_ply_GS.h"
#include <wrap/io_trimesh/export_ply.h>

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
                                           vcg::Use<MyEdge>     ::AsEdgeType,
                                           vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Color4b, vcg::vertex::Normal3f, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj,  vcg::face::Color4b, vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<   MyUsedTypes> {};

class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

using namespace std;


int main(int argc, char *argv[])
{
    if(argc != 10) {
        cout << "Insufficient arguments" << endl;
        cout << "Expected args: input.ply minXBox minYBox minZBox maxXBox maxYBox maxZBox theta phi" << endl;
        return -1;
    }

    MyMesh mEllips, mSampledEllips;
    bool isContained = true;

    // Setup box
    vcg::Point3<float> minBox = vcg::Point3<float>(stof(argv[2]), stof(argv[3]), stof(argv[4]));
    vcg::Point3<float> maxBox = vcg::Point3<float>(stof(argv[5]), stof(argv[6]), stof(argv[7]));
    vcg::Box3<float> box = vcg::Box3<float>(minBox, maxBox);
    // @TODO: could add support for box rotation

    int vIdx = 0;
    vcg::tri::io::PlyInfo pi;
    vcg::Matrix44f transf = vcg::Matrix44<float>();
    vcg::Quaternion<float> rotQuat;
    vcg::Point3<float> ellipseScale;
    vcg::Color4b ellipseColor;
    float theta = stof(argv[8]);
    float phi = stof(argv[9]);

    MyMesh gauss;
    const int DegreeSH = 1; // From 0 to 2
    MyMesh::PerVertexAttributeHandle<GaussianSplat<float,DegreeSH>> handleGauss = vcg::tri::io::ImporterPLYGS<MyMesh, DegreeSH>::Open(gauss, argv[1], pi);

    if(!vcg::tri::Allocator<MyMesh>::IsValidHandle(gauss, handleGauss))
    {
        return -1;
    }

    cout << "Number of Gaussian Splats: " << gauss.VN() << endl;

    // Loop through vertices, where each one represents a gaussian
    for(MyMesh::VertexIterator gi=gauss.vert.begin();gi!=gauss.vert.end();++gi) {
        isContained = true;
        if(vIdx % 1000 == 0) {
            cout << "GS id: " << vIdx << std::endl;
            cout << "Number of Gaussian Splats left: " << gauss.VN() << endl;
        }

        // Create solid
        vcg::tri::Octahedron(mEllips);

        // Scale
        vcg::tri::UpdatePosition<MyMesh>::Scale(mEllips, handleGauss[gi].getScale());
        // Rotate
        vcg::QuaternionToMatrix<float,vcg::Matrix44<float>>(handleGauss[gi].getRotation(), transf);
        vcg::tri::UpdatePosition<MyMesh>::Matrix(mEllips, transf);
        transf.SetIdentity();
        // Translate
        vcg::tri::UpdatePosition<MyMesh>::Translate(mEllips, gi->P());

        // Color
        vcg::tri::UpdateColor<MyMesh>::PerFaceConstant(mEllips, handleGauss[gi].getColor(theta, phi));

        // Delete vertices outside box
        for(MyMesh::VertexIterator vi=mEllips.vert.begin();vi!=mEllips.vert.end();++vi) {
            if(!box.IsIn(vi->P())) {
                vcg::tri::Allocator<MyMesh>::DeleteVertex(gauss, *gi);
                isContained = false;
                break;
            }
        }

        if(isContained) {
            // Add new ellipsoid to cluster mesh, for visual confirmation
            gi->SetS();
            vcg::tri::Append<MyMesh, MyMesh>::Mesh(mSampledEllips, mEllips);
        }

        mEllips.Clear();
        vIdx++;
    }

    // Save all remaining gaussians
    vcg::tri::Allocator<MyMesh>::CompactVertexVector(gauss);
    std::string filteredFileName = std::string("filteredGS.ply");
    vcg::tri::io::ExporterPLY<MyMesh>::Save(gauss, filteredFileName.c_str(), true, pi);
    // Delete all but the sampled ones
    vcg::tri::UpdateSelection<MyMesh>::VertexInvert(gauss);
    for(MyMesh::VertexIterator gi=gauss.vert.begin();gi!=gauss.vert.end();++gi) {
        if(gi->IsS())
            vcg::tri::Allocator<MyMesh>::DeleteVertex(gauss, *gi);
    }
    vcg::tri::Allocator<MyMesh>::CompactVertexVector(gauss);
    std::string gaussFileName = std::string("gaussSamples.ply");
    vcg::tri::io::ExporterPLY<MyMesh>::Save(gauss, gaussFileName.c_str(), true, pi);
    // Save ellipsoids from sampled vertices
    std::string ellipsFileName = std::string("ellipseSamplesSH") + std::to_string(DegreeSH) + ".ply";
    vcg::tri::io::ExporterPLY<MyMesh>::Save(mSampledEllips,ellipsFileName.c_str(), vcg::tri::io::Mask::IOM_FACECOLOR+vcg::tri::io::Mask::IOM_VERTCOLOR);
    mSampledEllips.Clear();

    return 0;
}
