#include <iostream>
#include <cmath>

#include <vcg/math/base.h>
#include <vcg/math/quaternion.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/space/fitting3.h>
#include <vcg/space/polygon3.h>
#include <vcg/space/point.h>
#include <vcg/math/matrix44.h>
#include <vcg/math/random_generator.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

#include <Eigen/src/Geometry/OrthoMethods.h>

#include "../trimesh_gaussian/export_ply_GS.h"

using namespace std;
using namespace vcg;

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>   ::AsVertexType,
                                      Use<MyEdge>     ::AsEdgeType,
                                      Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Color4b, vertex::Normal3f, vertex::VFAdj, vertex::BitFlags  >{};
class MyFace    : public Face<   MyUsedTypes, face::FFAdj, face::VFAdj, face::Normal3f, face::Color4b, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<   MyUsedTypes> {};

class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge> > {};

template <typename MeshType>
class GaussianSplatConverter {
public:
    
    typedef typename MeshType::VertexType     VertexType;
    typedef typename MeshType::VertexType::CoordType     CoordType;
    typedef typename MeshType::VertexPointer  VertexPointer;
    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename MeshType::FaceIterator FaceIterator;
    typedef typename MeshType::ScalarType			ScalarType;
    
    
 /**
 * @brief PerFaceSplatApprox
 * @param m the mesh to be converted into splats
 * @param m_gs the mesh to store the splats
 * 
 * This function converts a mesh into a set of gaussian splats, one for each face
 * placed at the barycenter of the face, oriented along the face normal, and scaled
 * by the size of the face.
 */
static void PerFaceSplat(MeshType &m, MeshType &m_gs)
{
    m_gs.Clear();
    tri::UpdateNormal<MeshType>::PerFaceNormalized(m);
    tri::UpdateColor<MeshType>::PerFaceFromVertex(m);
    auto handleGS = tri::Allocator<MeshType>::template AddPerVertexAttribute<GaussianSplat<float,1>>(m_gs, "gs");
    tri::Allocator<MeshType>::AddVertices(m_gs, m.face.size());
    
    VertexIterator vi = m_gs.vert.begin();
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi,++vi) {
        Point3f normal = fi->N();
        Point3f barycenter = Barycenter(*fi);
        Point3f v1 = fi->V(0)->cP() - barycenter;
        float radius = v1.Norm();
        Point3f v2 = normal^v1;
        
        v1.Normalize();
        v2.Normalize();
        vcg::Matrix44<float> rotMat;
        rotMat.SetColumn(0, normal);
        rotMat.SetColumn(0, v1);
        rotMat.SetColumn(0, v2);
        Quaternion<float> rotQuat;
        rotQuat.FromMatrix(rotMat);
        Quaternion<float> rotQuat2(1,0,0,0); // Rotation irrelevant for a sphere
        
        Point3f scale(radius/10, radius, radius);
        vi->P()=barycenter;
        handleGS[vi] = GaussianSplat<float,1>(rotQuat2, scale, fi->C() );
    }
}

};




template<typename _MeshType>
void computeVertexAvgDist(_MeshType& mesh, int nNeighbors)
{
    typedef typename _MeshType::ScalarType    ScalarType;
    if (!vcg::tri::HasPerVertexAttribute(mesh, "radius")) {
        vcg::tri::Allocator<_MeshType>::template AddPerVertexAttribute<ScalarType>(mesh, "avgDist");
    }

    typename _MeshType::template PerVertexAttributeHandle<ScalarType> h;
    h = vcg::tri::Allocator<_MeshType>::template FindPerVertexAttribute<ScalarType>(mesh, "avgDist");
    assert(vcg::tri::Allocator<_MeshType>::template IsValidHandle<ScalarType>(mesh, h));
    
    VertexConstDataWrapper<_MeshType> DW(mesh);
    KdTree<ScalarType> knn(DW);
    
    typename vcg::KdTree<ScalarType>::PriorityQueue pq;
    for (size_t i = 0; i < mesh.vert.size(); i++) {
        knn.doQueryK(mesh.vert[i].cP(), nNeighbors, pq);
        ScalarType totDist = 0;
        for(int j=0;j<pq.getNofElements();j++)
            totDist += sqrt(pq.getWeight(j));
        h[i] = totDist / ScalarType(pq.getNofElements());
    }
}

template<class _MeshType>
void computePCA(const vector<typename _MeshType::CoordType> &pointVec,
             typename _MeshType::CoordType PCA[])
{
    typedef typename _MeshType::CoordType CoordType;
    typedef typename _MeshType::ScalarType ScalarType;

    //compute the covariance matrix
    Eigen::Matrix<ScalarType,3,3> EigenCovMat;
    CoordType Barycenter;

    vcg::ComputeCovarianceMatrix<ScalarType>(pointVec, Barycenter, EigenCovMat);

    // Compute pca vectors
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ScalarType,3,3>> eig(EigenCovMat);

    Eigen::Matrix<ScalarType,1,3> eval = eig.eigenvalues();
    Eigen::Matrix<ScalarType,3,3> evec = eig.eigenvectors();

    eval = eval.cwiseAbs();
    int normInd,maxInd,minInd;

    ///get min and max coff ..
    ///the minumum is the Normal
    ///the other two the anisotropy directions
    eval.minCoeff(&normInd);
    eval.maxCoeff(&maxInd);
    minInd=(maxInd+1)%3;

    if (minInd==normInd)minInd=(normInd+1)%3;
    assert((minInd!=normInd)&&(minInd!=maxInd)&&(minInd!=maxInd));

    ///maximum direction of PCA
    PCA[0][0] = evec(0,maxInd);
    PCA[0][1] = evec(1,maxInd);
    PCA[0][2] = evec(2,maxInd);
    ///minimum direction of PCA
    PCA[1][0] = evec(0,minInd);
    PCA[1][1] = evec(1,minInd);
    PCA[1][2] = evec(2,minInd);
    ///Normal direction
    PCA[2][0] = evec(0,normInd);
    PCA[2][1] = evec(1,normInd);
    PCA[2][2] = evec(2,normInd);

    ScalarType LX=sqrt(eval[maxInd]);
    ScalarType LY=sqrt(eval[minInd]);
    //ScalarType LZ=sqrt(eval[normInd]);

    ///scale the directions
    PCA[0]*=LX;
    PCA[1]*=LY;
    //PCA[2]*=LZ;//.Normalize();
    PCA[2].Normalize();
}

template<typename _MeshType>
void computeVertexPCA(_MeshType& mesh, int nNeighbors)
{
    typedef typename _MeshType::ScalarType    Scalar;
    if (!vcg::tri::HasPerVertexAttribute(mesh, "pca")) {
        vcg::tri::Allocator<_MeshType>::template AddPerVertexAttribute<vector<vcg::Point3<Scalar>>>(mesh, "pca");
    }

    typename _MeshType::template PerVertexAttributeHandle<vector<vcg::Point3<Scalar>>> h;
    h = vcg::tri::Allocator<_MeshType>::template FindPerVertexAttribute<vector<vcg::Point3<Scalar>>>(mesh, "pca");
    assert(vcg::tri::Allocator<_MeshType>::template IsValidHandle<vector<vcg::Point3<Scalar>>>(mesh, h));

    auto positions = vcg::ConstDataWrapper<vcg::Point3<Scalar>>(
        &mesh.vert[0].P(),
        mesh.vert.size(),
        size_t(mesh.vert[1].P().V()) - size_t(mesh.vert[0].P().V()));

    vcg::KdTree<Scalar> knn(positions);
    typename vcg::KdTree<Scalar>::PriorityQueue pq;
    for (size_t i = 0; i < mesh.vert.size(); i++) {
        knn.doQueryK(mesh.vert[i].cP(), nNeighbors, pq);
        vector<Point3<Scalar>> pointsVec;
        for(int j=0;j<pq.getNofElements();j++)
        {
            pointsVec.push_back(mesh.vert[pq.getIndex(j)].cP());
        }
        h[i].resize(3);
        computePCA<_MeshType>(pointsVec, &h[i][0]);
    }
}
// Uniform splat approximation
// we generate a splat for each vertex
// all splats are spherical scaled by the average distance to the n nearest neighbours

void uniformSplatApprox(MyMesh &m)
{
    // Compute radius attribute of each vertex using the 10 nearest neighbours
    computeVertexAvgDist<MyMesh>(m, 6);

    MyMesh::PerVertexAttributeHandle<float> radiusH = tri::Allocator<MyMesh>::GetPerVertexAttribute<float>(m, "avgDist");
    MyMesh::PerVertexAttributeHandle<GaussianSplat<float,1>> handleGS = tri::Allocator<MyMesh>::AddPerVertexAttribute<GaussianSplat<float,1>>(m, "gs");

    for(MyMesh::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) {
        Quaternion<float> rotQuat(1,0,0,0); // Rotation irrelevant for a sphere
        Point3f scale(radiusH[vi], radiusH[vi], radiusH[vi]);
        handleGS[vi] = GaussianSplat<float,1>(rotQuat, scale, vi->cC());
    }
}

void flatSplatApprox(MyMesh &m)
{
    // Setup attributes
    MyMesh::PerVertexAttributeHandle<GaussianSplat<float,1>> handleGS = tri::Allocator<MyMesh>::AddPerVertexAttribute<GaussianSplat<float,1>>(m, "gs");

    // Compute normals per vertex
    vcg::tri::UpdateNormal<MyMesh>::PerVertex(m);

    for(MyMesh::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) {

        Point3f normal = vi->cN();
        vcg::Matrix44<float> rotMat;
        vcg::Matrix44<float> translMat1, translMat2;

        // Translate vertex to origin and translate normal accordingly
        translMat1.SetTranslate(-vi->cP());
        translMat2.SetTranslate(vi->cP());

        // Compute rotation
        Point3f newNormal = translMat1 * normal;

        float angleZ = Angle(Point3<float>(0,0,1), newNormal);
        Point3f orthZ = Point3<float>(0,0,1) ^ newNormal;
        vcg::Matrix44<float> matZ = rotMat.SetRotateRad(angleZ, orthZ);

        float size = 0.01;
        Quaternion<float> rotQuat;
        rotQuat.FromMatrix(translMat2*matZ*translMat1);
        Point3f scale(1, 1, 0.1);
        handleGS[vi] = GaussianSplat<float,1>(rotQuat,size*scale, vi->cC());
    }
}

void flatSplatApproxPCA(MyMesh &m)
{
    // Setup attributes
    MyMesh::PerVertexAttributeHandle<vector<Point3f>> handlePCA = tri::Allocator<MyMesh>::GetPerVertexAttribute<vector<Point3f>>(m, "pca");
    MyMesh::PerVertexAttributeHandle<GaussianSplat<float,3>> handleGS = tri::Allocator<MyMesh>::AddPerVertexAttribute<GaussianSplat<float,3>>(m, "gs");

    // Compute PCA per vertex
    computeVertexPCA<MyMesh>(m, 10);
    cout << "PCA computed" << endl;

    for(MyMesh::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) {
        vector<Point3f> pca = handlePCA[vi];
        vcg::Matrix33<float> rotMat;

        // Compute scale
        Point3f scale(pca[0].Norm(), pca[1].Norm(), pca[0].Norm()/20);

        // Compute rotation
        // Translate vertex point to origin, together with normal
        Point3f normal = vi->cN();
        Point3f newNormal = normal - vi->cP();

        // Align normal of vertex with normal of pca vectors
        float angleNormal = vcg::Angle(pca[2], newNormal);
        Point3f orthNormal = pca[2] ^ newNormal.normalized();
        vcg::Matrix33<float> matNormal = rotMat.SetRotateRad(angleNormal, orthNormal);

        // Align one of the pca axis by rotating around normal

        // Get angles between std axis and pca axis
        // Construct orthogonal vector to plane on which angle lies and rotate around it
        float angleMax = vcg::AngleN(pca[0].normalized(), Point3f(0,0,1));
        Point3f orthMax = Point3<float>(0,1,0) ^ pca[0].normalized();
        vcg::Matrix33<float> matMax = rotMat.SetRotateRad(angleMax, orthMax);

        float angleMin = vcg::AngleN(pca[1].normalized(), Point3f(0,1,0));
        Point3f orthMin = (Point3f(0,1,0)) ^ pca[1].normalized();
        vcg::Matrix33<float> matMin = rotMat.SetRotateRad(angleMin, orthMin);

        Quaternion<float> rotQuat;
        rotQuat.FromMatrix(matMin*matMax);
        handleGS[vi] = GaussianSplat<float,3>(rotQuat, scale, vi->cC());
    }
}

int main(int argc, char *argv[])
{
    if(argc != 3) {
        cout << "Insufficient arguments" << endl;
        cout << "Expected args: inputPointCloud.ply output.ply" << endl;
        return -1;
    }

    MyMesh inputMesh;
    MyMesh mPointCloud;

    // Load input mesh
    tri::io::Importer<MyMesh>::Open(inputMesh, argv[1]);

    //uniformSplatApprox(mPointCloud);
    // uniformSplatApprox(mPointCloud);
    GaussianSplatConverter<MyMesh>::PerFaceSplat(inputMesh, mPointCloud);
    //flatSplatApproxPCA(mPointCloud);

    tri::io::ExporterPLYGS<MyMesh, 1>::Save(mPointCloud, argv[2], true);

    return 0;
}
