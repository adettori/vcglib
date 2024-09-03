#include <iostream>
#include <algorithm>
#include <cmath>
#include <vcg/math/base.h>
#include <vcg/math/quaternion.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/intersection.h>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
                                           vcg::Use<MyEdge>     ::AsEdgeType,
                                           vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Color4b, vcg::vertex::Normal3f, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj,  vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<   MyUsedTypes> {};

class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

const std::string properties[] = {"nxx", "ny", "nz", "f_dc_0", "f_dc_1", "f_dc_2", "f_rest_0", "f_rest_1", "f_rest_2", "f_rest_3", "f_rest_4", "f_rest_5", "f_rest_6", "f_rest_7", "f_rest_8", "f_rest_9", "f_rest_10", "f_rest_11", "f_rest_12", "f_rest_13", "f_rest_14", "f_rest_15", "f_rest_16", "f_rest_17", "f_rest_18", "f_rest_19", "f_rest_20", "f_rest_21", "f_rest_22", "f_rest_23", "f_rest_24", "f_rest_25", "f_rest_26", "f_rest_27", "f_rest_28", "f_rest_29", "f_rest_30", "f_rest_31", "f_rest_32", "f_rest_33", "f_rest_34", "f_rest_35", "f_rest_36", "f_rest_37", "f_rest_38", "f_rest_39", "f_rest_40", "f_rest_41", "f_rest_42", "f_rest_43", "f_rest_44", "opacity", "scale_0", "scale_1", "scale_2", "rot_0", "rot_1", "rot_2", "rot_3"};

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 9) {
        cout << "Insufficient arguments" << endl;
        cout << "Expected args: input.ply output.ply minXBox minYBox minZBox maxXBox maxYBox maxZBox" << endl;
        return -1;
    }

    int propLen = sizeof(properties)/sizeof(properties[0]);
    MyMesh mGauss, mEllips, mEllipsCluster, mBox;
    bool isContained = true;
    vcg::tri::io::PlyInfo piOpen;
    MyMesh::PerVertexAttributeHandle<float>* handler_arr = new MyMesh::PerVertexAttributeHandle<float>[sizeof(properties)];

    for(int i=0;i<propLen;i++) {
        piOpen.AddPerVertexFloatAttribute(properties[i]);
        // add a per-vertex attribute with type float
        handler_arr[i] = vcg::tri::Allocator<MyMesh>:: GetPerVertexAttribute<float> (mGauss, properties[i]);
    }

    // Setup box
    vcg::Point3<float> minBox = vcg::Point3<float>(stof(argv[3]), stof(argv[4]), stof(argv[5]));
    vcg::Point3<float> maxBox = vcg::Point3<float>(stof(argv[6]), stof(argv[7]), stof(argv[8]));
    vcg::Box3<float> box = vcg::Box3<float>(minBox, maxBox);
    vcg::tri::Box<MyMesh>(mBox, box);
    // @TODO: could add support for box rotation
    vcg::tri::io::ExporterPLY<MyMesh>::Save(mBox,"box.ply");

    // Find indices of relevant properties
    // Scale
    int propScaleX = find(&properties[0], properties + propLen, "scale_0") - properties;
    int propScaleY = find(&properties[0], properties + propLen, "scale_1") - properties;
    int propScaleZ = find(&properties[0], properties + propLen, "scale_2") - properties;
    // Rotation
    int propRot0 = find(&properties[0], properties + propLen, "rot_0") - properties;
    int propRot1 = find(&properties[0], properties + propLen, "rot_1") - properties;
    int propRot2 = find(&properties[0], properties + propLen, "rot_2") - properties;
    int propRot3 = find(&properties[0], properties + propLen, "rot_3") - properties;

    int ret = vcg::tri::io::ImporterPLY<MyMesh>::Open(mGauss, argv[1], piOpen);
    if(ret!=0)
    {
        printf("Unable to open %s for '%s'\n", argv[1], vcg::tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
        return -1;
    }

    cout << "Number of Gaussian Splats: " << mGauss.VN() << endl;

    int vIdx = 0;
    vcg::Matrix44f transf = vcg::Matrix44<float>();
    vcg::Quaternion<float> rotQuat;
    vcg::Point3<float> ellipsePos;
    vcg::Point3<float> ellipseScale;

    // Loop through vertices, where each one represents a gaussian
    for(MyMesh::VertexIterator gi=mGauss.vert.begin();gi!=mGauss.vert.end();++gi) {
        isContained = true;
        if(vIdx % 1000 == 0) {
            cout << "GS id: " << vIdx << std::endl;
            cout << "Number of Gaussian Splats left: " << mGauss.VN() << endl;
        }

        // Create solid
        vcg::tri::SuperEllipsoid(mEllips, 2, 2, 2, 24, 12); // r=s=t=2 to get an ellipsoid

        // Rotate
        rotQuat = vcg::Quaternion<float>(handler_arr[propRot0][gi], handler_arr[propRot1][gi], handler_arr[propRot2][gi], handler_arr[propRot3][gi]);
        rotQuat.ToMatrix(transf);
        vcg::tri::UpdatePosition<MyMesh>::Matrix(mEllips, transf);
        transf.SetIdentity();

        // Translate and scale
        ellipsePos = vcg::Point3<float>(gi->P().X(), gi->P().Y(), gi->P().Z());
        // Need to exponentiate the scale values
        ellipseScale = vcg::Point3<float>(exp(handler_arr[propScaleX][gi]), exp(handler_arr[propScaleY][gi]), exp(handler_arr[propScaleZ][gi]));
        transf.SetTranslate(ellipsePos);
        transf.SetScale(ellipseScale);
        vcg::tri::UpdatePosition<MyMesh>::Matrix(mEllips, transf);

        // Delete vertices outside box
        for(MyMesh::VertexIterator vi=mEllips.vert.begin();vi!=mEllips.vert.end();++vi) {
            if(!box.IsIn(vi->P())) {
                //cout << box.min.X() << " " << box.min.Y() << " " << box.min.Z() << endl;
                //cout << vi->P().X() << " " << vi->P().Y() << " " << vi->P().Z() << endl; // Importer not getting the correct attributes from file
                //cout << box.max.X() << " " << box.max.Y() << " " << box.max.Z() << endl;
                vcg::tri::Allocator<MyMesh>::DeleteVertex(mGauss, *gi);
                isContained = false;
                break;
            }
        }

        if(vIdx % 1000 == 0 && isContained)
            // Add new ellipsoid to cluster mesh, for visual confirmation
            vcg::tri::Append<MyMesh, MyMesh>::Mesh(mEllipsCluster, mEllips);

        mEllips.Clear();
        vIdx++;
    }

    vcg::tri::Allocator<MyMesh>::CompactVertexVector(mGauss);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(mGauss,"filteredGs.ply", true, piOpen);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(mEllipsCluster,"gaussCluster.ply");
    mEllipsCluster.Clear();
    delete[] handler_arr;

    return 0;
}
