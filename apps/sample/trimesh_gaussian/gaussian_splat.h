#ifndef GAUSSIAN_SPLAT_H
#define GAUSSIAN_SPLAT_H

#include <cmath>
#include <set>
#include <vcg/math/base.h>
#include <vcg/math/quaternion.h>
#include <wrap/io_trimesh/io_ply.h>

using namespace std;

template <typename MeshType, typename ScalarType> class GaussianSplat {

public:
    typedef typename MeshType::VertexType     VertexType;
    typedef typename MeshType::VertexPointer  VertexPointer;
    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename MeshType::VertContainer VertContainer;

    typedef typename MeshType::EdgeType     EdgeType;
    typedef typename MeshType::EdgePointer  EdgePointer;
    typedef typename MeshType::EdgeIterator EdgeIterator;
    typedef typename MeshType::EdgeContainer EdgeContainer;

    typedef typename MeshType::FaceType       FaceType;
    typedef typename MeshType::FacePointer    FacePointer;
    typedef typename MeshType::FaceIterator   FaceIterator;
    typedef typename MeshType::FaceContainer FaceContainer;

    typedef typename MeshType::HEdgeType     HEdgeType;
    typedef typename MeshType::HEdgePointer  HEdgePointer;
    typedef typename MeshType::HEdgeIterator HEdgeIterator;
    typedef typename MeshType::HEdgeContainer HEdgeContainer;

    typedef typename MeshType::TetraType TetraType;
    typedef typename MeshType::TetraPointer TetraPointer;
    typedef typename MeshType::TetraIterator TetraIterator;
    typedef typename MeshType::TetraContainer TetraContainer;

    typedef typename MeshType::CoordType     CoordType;

    typedef typename MeshType::PointerToAttribute PointerToAttribute;
    typedef typename std::set<PointerToAttribute>::iterator AttrIterator;
    typedef typename std::set<PointerToAttribute>::const_iterator AttrConstIterator;
    typedef typename std::set<PointerToAttribute >::iterator PAIte;

    int degreeSH;
    vcg::Quaternion<ScalarType> rot;
    vcg::Point3<ScalarType> scale;

    // Constructors
    GaussianSplat() {}

    GaussianSplat(vcg::Quaternion<ScalarType> rot, vcg::Point3<ScalarType> scale, int degreeSH=0) {
        this->rot = rot;
        this->scale = scale;
        this->degreeSH = degreeSH;
    }

    ~GaussianSplat() {
    }
};

#endif // GAUSSIAN_SPLAT_H
