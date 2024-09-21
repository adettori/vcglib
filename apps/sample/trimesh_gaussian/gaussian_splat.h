#ifndef GAUSSIAN_SPLAT_H
#define GAUSSIAN_SPLAT_H

#include <cmath>
#include <vcg/math/base.h>
#include <vcg/math/quaternion.h>
#include <wrap/io_trimesh/io_ply.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/math/spherical_harmonics.h>

using namespace std;

template <typename ScalarType, int DegreeSH>
class GaussianSplat {

private:
    vcg::Color4b color; // SH0 color
    vcg::Quaternion<ScalarType> rot;
    vcg::Point3<ScalarType> scale;
    vcg::math::SphericalHarmonics<ScalarType, DegreeSH+1> sphR, sphG, sphB;

    static int clamp(int v, int lo, int hi) {
        return min(hi, max(lo, v));
    }

    static int fdcToColor(float value)
    {
        return clamp(round((0.5 + value)*255), 0, 255);
    }

public:
    constexpr static const float SH_C0 = 0.28209479177387814;

    // Constructors
    GaussianSplat() {}

    GaussianSplat(vcg::Quaternion<ScalarType> rot, vcg::Point3<ScalarType> scale,
                  std::vector<ScalarType> &vecSHR, std::vector<ScalarType> &vecSHG, std::vector<ScalarType> &vecSHB,ScalarType alpha) {
        this->rot = rot;
        this->scale = scale;

        this->sphR = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&vecSHR[0]);  // red
        this->sphG = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&vecSHG[0]);  // green
        this->sphB = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&vecSHB[0]);  // blue

        // Set color
        this->color = vcg::Color4b(
            fdcToColor(SH_C0*vecSHR[0]),
            fdcToColor(SH_C0*vecSHG[0]),
            fdcToColor(SH_C0*vecSHB[0]),
            clamp(round(1 / (1 + exp(-alpha)) * 255), 0, 255)
        );
    }

    GaussianSplat(vcg::Quaternion<ScalarType> rot, vcg::Point3<ScalarType> scale, vcg::Color4b color) {
        this->rot = rot;
        this->scale = scale;
        this->color = color;
    }

    ~GaussianSplat() {
    }

    vcg::Color4b getColor()
    {
        return this->color;
    }

    vcg::Color4b getColor(ScalarType theta, ScalarType phi)
    {
        return vcg::Color4b(fdcToColor(sphR(theta, phi)), // r
                            fdcToColor(sphG(theta, phi)), // g
                            fdcToColor(sphB(theta, phi)), // b
                            this->colorValues[3]); // alpha
    }

    vcg::Quaternion<ScalarType> getRotation()
    {
        return this->rot;
    }

    vcg::Point3<ScalarType> getScale()
    {
        return this->scale;
    }
};

#endif // GAUSSIAN_SPLAT_H
