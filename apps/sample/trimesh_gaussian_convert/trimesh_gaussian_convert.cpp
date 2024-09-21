#include <iostream>
#include <cmath>
#include <vcg/math/base.h>
#include <vcg/math/quaternion.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>
#include "../trimesh_gaussian/import_ply_GS.h"
using namespace std;
using namespace vcg;

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>   ::AsVertexType,
                                      Use<MyEdge>     ::AsEdgeType,
                                      Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Color4b, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face<   MyUsedTypes, face::FFAdj,  face::Color4b, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<   MyUsedTypes> {};

class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge> > {};

template<typename _MeshType>
void computeVertexRadius(_MeshType& mesh, int nNeighbors)
{
    typedef typename _MeshType::ScalarType    Scalar;
    if (!vcg::tri::HasPerVertexAttribute(mesh, "radius")) {
        vcg::tri::Allocator<_MeshType>::template AddPerVertexAttribute<Scalar>(mesh, "radius");
    }

    typename _MeshType::template PerVertexAttributeHandle<Scalar> h;
    h = vcg::tri::Allocator<_MeshType>::template FindPerVertexAttribute<Scalar>(mesh, "radius");
    assert(vcg::tri::Allocator<_MeshType>::template IsValidHandle<Scalar>(mesh, h));

    auto positions = vcg::ConstDataWrapper<vcg::Point3<Scalar>>(
        &mesh.vert[0].P(),
        mesh.vert.size(),
        size_t(mesh.vert[1].P().V()) - size_t(mesh.vert[0].P().V()));

    vcg::KdTree<Scalar> knn(positions);
    typename vcg::KdTree<Scalar>::PriorityQueue pq;
    for (size_t i = 0; i < mesh.vert.size(); i++) {
        knn.doQueryK(mesh.vert[i].cP(), nNeighbors, pq);
        h[i] = 2. * sqrt(pq.getTopWeight() / float(pq.getNofElements()));
    }
}

int main(int argc, char *argv[])
{
    if(argc < 3) {
        cout << "Insufficient arguments" << endl;
        cout << "Expected args: input.ply output.ply" << endl;
        return -1;
    }

    MyMesh m;

    // Load input mesh
    tri::io::Importer<MyMesh>::Open(m, argv[1]);

    // Compute radius attribute of each vertex using the 10 nearest neighbours
    computeVertexRadius<MyMesh>(m, 10);

    MyMesh::PerVertexAttributeHandle<float> handleR = tri::Allocator<MyMesh>::GetPerVertexAttribute<float>(m, "radius");
    MyMesh::PerVertexAttributeHandle<GaussianSplat<float,0>> handleGS = tri::Allocator<MyMesh>::AddPerVertexAttribute<GaussianSplat<float,0>>(m, "gs");
    vector<float> shR, shG, shB;

    for(MyMesh::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) {
        Quaternion<float> rotQuat;
        float scaleValue = handleR[vi];
        vcg::Point3<float> scale(scaleValue, scaleValue, scaleValue);
        handleGS[vi] = GaussianSplat<float,0>(rotQuat, scale, vi->cC());
    }

    return 0;
}
