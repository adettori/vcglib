#ifndef GAUSSIAN_SPLAT_H
#define GAUSSIAN_SPLAT_H

#include <cmath>
#include <vcg/math/base.h>
#include <vcg/math/quaternion.h>
#include <wrap/io_trimesh/io_ply.h>
#include <vcg/math/spherical_harmonics.h>

using namespace std;

template <typename ScalarType> class GaussianSplat {

public:

    int degreeSH;
    vcg::Quaternion<ScalarType> rot;
    vcg::Point3<ScalarType> scale;
    ScalarType listSH[];

    // Constructors
    GaussianSplat() {}

    GaussianSplat(vcg::Quaternion<ScalarType> rot, vcg::Point3<ScalarType> scale, ScalarType listSH[], int degreeSH=0) {
        this->rot = rot;
        this->scale = scale;
        this->degreeSH = degreeSH;
        std::copy(listSH, listSH, this->listSH);
    }

    ~GaussianSplat() {
    }
};

#endif // GAUSSIAN_SPLAT_H
