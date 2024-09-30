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

private:

static void computeVertexAvgDist(MeshType& mesh, int nNeighbors)
{
    typedef typename MeshType::ScalarType    ScalarType;
    if (!vcg::tri::HasPerVertexAttribute(mesh, "radius")) {
        vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<ScalarType>(mesh, "avgDist");
    }

    typename MeshType::template PerVertexAttributeHandle<ScalarType> h;
    h = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<ScalarType>(mesh, "avgDist");
    assert(vcg::tri::Allocator<MeshType>::template IsValidHandle<ScalarType>(mesh, h));

    VertexConstDataWrapper<MeshType> DW(mesh);
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

static void computePCA(const vector<typename MeshType::CoordType> &pointVec,
                typename MeshType::CoordType PCA[])
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

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
    PCA[2].Normalize();
}

static void computeVertexPCA(MeshType& mesh, int nNeighbors)
{
    typedef typename MeshType::ScalarType    Scalar;
    if (!vcg::tri::HasPerVertexAttribute(mesh, "pca")) {
        vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<vector<vcg::Point3<Scalar>>>(mesh, "pca");
    }

    typename MeshType::template PerVertexAttributeHandle<vector<vcg::Point3<Scalar>>> h;
    h = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<vector<vcg::Point3<Scalar>>>(mesh, "pca");
    assert(vcg::tri::Allocator<MeshType>::template IsValidHandle<vector<vcg::Point3<Scalar>>>(mesh, h));

    VertexConstDataWrapper<MeshType> DW(mesh);

    vcg::KdTree<Scalar> knn(DW);
    typename vcg::KdTree<Scalar>::PriorityQueue pq;
    for (size_t i = 0; i < mesh.vert.size(); i++) {
        knn.doQueryK(mesh.vert[i].cP(), nNeighbors, pq);
        vector<Point3<Scalar>> pointsVec;
        for(int j=0;j<pq.getNofElements();j++)
        {
            pointsVec.push_back(mesh.vert[pq.getIndex(j)].cP());
        }
        h[i].resize(3);
        computePCA(pointsVec, &h[i][0]);
    }
}

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

// Uniform splat approximation
// we generate a splat for each vertex
// all splats are spherical scaled by the average distance to the n nearest neighbours

static void uniformSplatApprox(MyMesh &m)
{
    // Compute radius attribute of each vertex using the 10 nearest neighbours
    computeVertexAvgDist(m, 6);

    MyMesh::PerVertexAttributeHandle<float> radiusH = tri::Allocator<MyMesh>::GetPerVertexAttribute<float>(m, "avgDist");
    MyMesh::PerVertexAttributeHandle<GaussianSplat<float,1>> handleGS = tri::Allocator<MyMesh>::AddPerVertexAttribute<GaussianSplat<float,1>>(m, "gs");

    for(MyMesh::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) {
        Quaternion<float> rotQuat(1,0,0,0); // Rotation irrelevant for a sphere
        Point3f scale(radiusH[vi], radiusH[vi], radiusH[vi]);
        handleGS[vi] = GaussianSplat<float,1>(rotQuat, scale, vi->cC());
    }
}

static void flatSplatApprox(MyMesh &m)
{
    // Setup attributes
    MyMesh::PerVertexAttributeHandle<GaussianSplat<float,1>> handleGS = tri::Allocator<MyMesh>::AddPerVertexAttribute<GaussianSplat<float,1>>(m, "gs");

    // Compute normals per vertex
    vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);

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

static void flatSplatApproxPCA(MyMesh &m)
{
    // Setup attributes
    MyMesh::PerVertexAttributeHandle<vector<Point3f>> handlePCA = tri::Allocator<MyMesh>::GetPerVertexAttribute<vector<Point3f>>(m, "pca");
    MyMesh::PerVertexAttributeHandle<GaussianSplat<float,1>> handleGS = tri::Allocator<MyMesh>::AddPerVertexAttribute<GaussianSplat<float,1>>(m, "gs");

    // Compute PCA per vertex
    computeVertexPCA(m, 5);
    cout << "PCA computed" << endl;

    for(MyMesh::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) {
        vector<Point3f> pca = handlePCA[vi];
        vcg::Matrix44<float> rotMat;

        // Compute scale
        Point3f scale(pca[0].Norm(), pca[1].Norm(), pca[0].Norm()/20);

        // Compute rotation
        // Translate vertex point to origin, together with normal
        Point3f normal = vi->cN();

        // Build orthogonal basis starting from newNormal
        Eigen::Vector3f normEigVec;
        Eigen::Vector3f orthEigVec;
        normal.ToEigenVector(normEigVec);
        orthEigVec = normEigVec.unitOrthogonal();
        Point3f normalBasis1;
        normalBasis1.FromEigenVector(orthEigVec);
        Point3f normalBasis2 = rotMat.SetRotateDeg(-90, normal) * normalBasis1;

        // Align normal of vertex with normal of pca vectors
        float angleNormal = vcg::Angle(normal, pca[2]);
        Point3f orthNormal = pca[2] ^ normal;
        vcg::Matrix44<float> matNormal = rotMat.SetRotateRad(-angleNormal, orthNormal);

        // Update normal basis
        normal = matNormal * normal;
        normalBasis1 = matNormal * normalBasis1;
        normalBasis2 = matNormal * normalBasis2;

        // Align one of the pca axis by rotating around normal
        float angleMax = vcg::Angle(normalBasis1, pca[0]);
        vcg::Matrix44<float> matMax = rotMat.SetRotateRad(-angleMax, normal);

        // Update normal basis
        normalBasis1 = matMax * normalBasis1;
        normalBasis2 = matMax * normalBasis2;

        float angleMin = vcg::Angle(normalBasis2, pca[1]);
        vcg::Matrix44<float> matMin = rotMat.SetRotateRad(-angleMin, normal);

        /*
        // Update normal basis (not used, just to make procedure clear)
        normalBasis1 = matMin * normalBasis1;
        normalBasis2 = matMin * normalBasis2;

        // Debug info
        cout << normal[0] << " " << normal[1] << " " << normal[2] << endl;
        cout << normalBasis1[0] << " " << normalBasis1[1] << " " << normalBasis1[2] << endl;
        cout << normalBasis2[0] << " " << normalBasis2[1] << " " << normalBasis2[2] << endl;
        cout << pca[2][0] << " " << pca[2][1] << " " << pca[2][2] << endl;
        cout << normalBasis1 * normalBasis2 << endl;
        cout << vcg::Angle(normal, normalBasis1) << endl;
        cout << vcg::Angle(normal, normalBasis2) << endl;
        cout << vcg::Angle(normal, pca[0]) << endl;
        cout << vcg::Angle(normal, pca[1]) << endl;
        cout << vcg::Angle(normal, pca[2]) << endl;
        cout << vcg::Angle(normalBasis1, pca[0]) << endl;
        cout << vcg::Angle(normalBasis2, pca[1]) << endl;
        cout << endl;
        exit(0);
        */

        Quaternion<float> rotQuat;
        rotQuat.FromMatrix(matMin*matMax*matNormal);
        handleGS[vi] = GaussianSplat<float,1>(rotQuat, scale, vi->cC());
    }
}

};

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
    tri::io::Importer<MyMesh>::Open(mPointCloud, argv[1]);

    // GaussianSplatConverter<MyMesh>::uniformSplatApprox(mPointCloud);
    // GaussianSplatConverter<MyMesh>::PerFaceSplat(inputMesh, mPointCloud);
    GaussianSplatConverter<MyMesh>::flatSplatApproxPCA(mPointCloud);

    tri::io::ExporterPLYGS<MyMesh, 1>::Save(mPointCloud, argv[2], true);

    return 0;
}
