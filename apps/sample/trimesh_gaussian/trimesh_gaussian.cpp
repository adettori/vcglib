#include <iostream>
#include <cmath>
#include <vcg/math/base.h>
#include <vcg/math/quaternion.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>
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
    if(argc != 9) {
        cout << "Insufficient arguments" << endl;
        cout << "Expected args: input.ply output.ply minXBox minYBox minZBox maxXBox maxYBox maxZBox" << endl;
        return -1;
    }

    MyMesh mEllips, mEllipsCluster, mBox;
    bool isContained = true;

    // Setup box
    vcg::Point3<float> minBox = vcg::Point3<float>(stof(argv[3]), stof(argv[4]), stof(argv[5]));
    vcg::Point3<float> maxBox = vcg::Point3<float>(stof(argv[6]), stof(argv[7]), stof(argv[8]));
    vcg::Box3<float> box = vcg::Box3<float>(minBox, maxBox);
    vcg::tri::Box<MyMesh>(mBox, box);
    // @TODO: could add support for box rotation
    vcg::tri::io::ExporterPLY<MyMesh>::Save(mBox,"box.ply");

    MyMesh gauss;
    MyMesh::PerVertexAttributeHandle<GaussianSplat<MyMesh, float>> handleGauss = vcg::tri::io::ImporterPLYGS<MyMesh, float>::Open(gauss, argv[1]);

    cout << "Number of Gaussian Splats: " << gauss.VN() << endl;

    int vIdx = 0;
    vcg::Matrix44f transf = vcg::Matrix44<float>();
    vcg::Quaternion<float> rotQuat;
    vcg::Point3<float> ellipsePos;
    vcg::Point3<float> ellipseScale;
    vcg::Color4b ellipseColor;

    // Loop through vertices, where each one represents a gaussian
    for(MyMesh::VertexIterator gi=gauss.vert.begin();gi!=gauss.vert.end();++gi) {
        isContained = true;
        if(vIdx % 1000 == 0) {
            cout << "GS id: " << vIdx << std::endl;
            cout << "Number of Gaussian Splats left: " << gauss.VN() << endl;
        }

        // Create solid
        vcg::tri::SuperEllipsoid(mEllips, 2, 2, 2, 24, 12); // r=s=t=2 to get an ellipsoid

        // Rotate
        handleGauss[gi].rot.ToMatrix(transf);
        vcg::tri::UpdatePosition<MyMesh>::Matrix(mEllips, transf);
        transf.SetIdentity();

        // Translate and scale
        ellipsePos = vcg::Point3<float>(gi->P().X(), gi->P().Y(), gi->P().Z());
        // Need to exponentiate the scale values
        ellipseScale = handleGauss[gi].scale;
        transf.SetTranslate(ellipsePos);
        transf.SetScale(ellipseScale);
        vcg::tri::UpdatePosition<MyMesh>::Matrix(mEllips, transf);

        // Color
        vcg::tri::UpdateColor<MyMesh>::PerFaceConstant(mEllips, gi->C());

        // Delete vertices outside box
        for(MyMesh::VertexIterator vi=mEllips.vert.begin();vi!=mEllips.vert.end();++vi) {
            if(!box.IsIn(vi->P())) {
                vcg::tri::Allocator<MyMesh>::DeleteVertex(gauss, *gi);
                isContained = false;
                break;
            }
        }

        if(vIdx % 1000 == 0 && isContained) {
            // Add new ellipsoid to cluster mesh, for visual confirmation
            //cout << int(mEllips.face[0].C().X()) << " " << int(mEllips.face[0].C().Y()) << " " << int(mEllips.face[0].C().Z()) << " " << int(mEllips.face[0].C().W()) << endl;
            vcg::tri::Append<MyMesh, MyMesh>::Mesh(mEllipsCluster, mEllips);
        }

        mEllips.Clear();
        vIdx++;
    }

    vcg::tri::Allocator<MyMesh>::CompactVertexVector(gauss);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(gauss, argv[2], true);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(mEllipsCluster,"gaussCluster.ply", vcg::tri::io::Mask::IOM_FACECOLOR+vcg::tri::io::Mask::IOM_VERTCOLOR);
    mEllipsCluster.Clear();

    return 0;
}
