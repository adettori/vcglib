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

template <typename MeshType, int DegreeSH>
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
    typedef typename MeshType::ScalarType Scalar;
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
    auto handleGS = tri::Allocator<MeshType>::template AddPerVertexAttribute<GaussianSplat<float,DegreeSH>>(m_gs, "gs");
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
        rotMat.SetColumn(1, v1);
        rotMat.SetColumn(2, v2);
        Quaternion<float> rotQuat;
        rotQuat.FromMatrix(rotMat);
        
        Point3f scale(radius/10, radius, radius);
        vi->P()=barycenter;
        handleGS[vi] = GaussianSplat<float,DegreeSH>(rotQuat, scale, fi->C());
    }
}

/**
 * @brief PerVertexUniformSplat
 * @param m the mesh to be converted into splats
 *
 * This function converts a mesh into a set of gaussian splats, one for each vertex,
 * with each splat being a sphere scaled according to the avg distance of its numNeigh neighbours.
 */
static void PerVertexUniformSplat(MeshType &m, int numNeigh)
{
    // Compute radius attribute of each vertex using the 10 nearest neighbours
    computeVertexAvgDist(m, numNeigh);

    auto radiusH = tri::Allocator<MeshType>::template GetPerVertexAttribute<float>(m, "avgDist");
    auto handleGS = tri::Allocator<MeshType>::template AddPerVertexAttribute<GaussianSplat<float,DegreeSH>>(m, "gs");

    for(typename MeshType::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) {
        Quaternion<float> rotQuat(1,0,0,0); // Rotation irrelevant for a sphere
        Point3f scale(radiusH[vi], radiusH[vi], radiusH[vi]);
        handleGS[vi] = GaussianSplat<float,DegreeSH>(rotQuat, scale, vi->C());
    }
}

/**
 * @brief PerVertexFlatSplat
 * @param m the mesh to be converted into splats
 *
 * This function converts a mesh into a set of gaussian splats, one for each vertex,
 * with each splat being a flattened sphere, scaled according to the input parameter size.
 */
static void PerVertexFlatSplat(MeshType &m, float size=0.01)
{
    // Setup attributes
    auto handleGS = tri::Allocator<MeshType>::template AddPerVertexAttribute<GaussianSplat<float,DegreeSH>>(m, "gs");

    for(typename MeshType::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) {

        Point3f normal = vi->cN();
        vcg::Matrix44<float> rotMat, rotMat1;

        // Compute rotation
        // Build orthogonal basis starting from normal
        Eigen::Vector3f normEigVec;
        Eigen::Vector3f orthEigVec;
        normal.ToEigenVector(normEigVec);
        orthEigVec = normEigVec.unitOrthogonal();
        Point3f normalBasis1;
        normalBasis1.FromEigenVector(orthEigVec);
        vcg::Matrix44<float> normalMat = rotMat.SetRotateDeg(90, normal);
        Point3f normalBasis2 = normalMat * normalBasis1;

        // Set rotation matrix
        rotMat.SetColumn(0, normalBasis1);
        rotMat.SetColumn(1, normalBasis2);
        rotMat.SetColumn(2, normal);

        Quaternion<float> rotQuat;
        rotQuat.FromMatrix(rotMat);
        Point3f scale(size, size, size/10);
        handleGS[vi] = GaussianSplat<float,DegreeSH>(rotQuat, scale, vi->C());
    }
}

/**
 * @brief PerVertexFlatSplatPCA
 * @param m the mesh to be converted into splats
 *
 * This function converts a mesh into a set of gaussian splats, one for each vertex,
 * with each splat being a flattened ellipsoid, scaled by the norms of the PCA vectors
 * computed on the numNeigh neighbours.
 */
static void PerVertexFlatSplatPCA(MeshType &m, int numNeigh=5)
{
    // Setup attributes
    auto handlePCA = tri::Allocator<MeshType>::template GetPerVertexAttribute<vector<Point3f>>(m, "pca");
    auto handleGS = tri::Allocator<MeshType>::template AddPerVertexAttribute<GaussianSplat<float,DegreeSH>>(m, "gs");

    // Compute PCA per vertex
    computeVertexPCA(m, numNeigh);
    cout << "PCA computed" << endl;

    for(typename MeshType::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) {
        vector<Point3f> pca = handlePCA[vi];
        vcg::Matrix44<float> rotMat;

        // Compute scale
        Point3f scale(pca[0].Norm(), pca[1].Norm(), pca[0].Norm()/20);

        // Build rotation matrix from axis (pca)
        rotMat.SetColumn(0, pca[0].normalized());
        rotMat.SetColumn(1, pca[1].normalized());
        rotMat.SetColumn(2, pca[2].normalized());

        Quaternion<float> rotQuat;
        rotQuat.FromMatrix(rotMat);
        handleGS[vi] = GaussianSplat<float,DegreeSH>(rotQuat, scale, vi->C());
    }
}

};

int main(int argc, char *argv[])
{
    if(argc != 3) {
        cout << "Insufficient arguments" << endl;
        cout << "Expected args: inputMesh.ply output.ply" << endl;
        return -1;
    }

    MyMesh inputMesh;
    MyMesh mPointCloud;
    const int DegreeSH = 0;

    // Load input mesh
    tri::io::Importer<MyMesh>::Open(inputMesh, argv[1]);

    // Chosen method
    //GaussianSplatConverter<MyMesh,DegreeSH>::PerFaceSplat(inputMesh, mPointCloud);
    //GaussianSplatConverter<MyMesh,DegreeSH>::PerVertexUniformSplat(inputMesh, 10);
    //GaussianSplatConverter<MyMesh,DegreeSH>::PerVertexFlatSplat(inputMesh);
    GaussianSplatConverter<MyMesh,DegreeSH>::PerVertexFlatSplatPCA(inputMesh);

    tri::io::ExporterPLYGS<MyMesh,DegreeSH>::Save(inputMesh, argv[2], true);

    return 0;
}
