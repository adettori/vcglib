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
    vcg::Color4b color; // Colors + opacity values
    vcg::Point4<ScalarType> colorValues; // SH0 + alpha values
    vcg::Quaternion<ScalarType> rot;
    vcg::Point3<ScalarType> scale;
    vcg::math::SphericalHarmonics<ScalarType, DegreeSH+1> sphR, sphG, sphB;

    static int clamp(int v, int lo, int hi) {
        return min(hi, max(lo, v));
    }

    static int fdcToColor(ScalarType value)
    {
        return clamp(round((0.5 + value)*255), 0, 255);
    }

    static ScalarType colorToFdc(ScalarType value)
    {
        return value/255-0.5;
    }

    static int opacityToAlpha(ScalarType opacity)
    {
        return clamp(round(1 / (1 + exp(-opacity)) * 255), 0, 255);
    }

    static int alphaToOpacity(int alpha)
    {
        ScalarType eps = 0.0000001;
        return log(alpha/(255-alpha+eps));
    }

    static vector<vector<ScalarType>> splitSHChannels(vector<ScalarType> &vecSH)
    {
        vector<vector<ScalarType>> channelSH = {vector<ScalarType>(), vector<ScalarType>(), vector<ScalarType>()}; // 1 vector per channel
        int elemsDegree, prevElems;
        // SH0
        channelSH[0].push_back(vecSH[0]);
        channelSH[1].push_back(vecSH[1]);
        channelSH[2].push_back(vecSH[2]);

        if(DegreeSH >= 1){
            prevElems = 3;
            elemsDegree = 3;
            for(int i=0;i<elemsDegree;i++)
            {
                for(int channel=0;channel<3;channel++)
                    channelSH[channel].push_back(vecSH[prevElems+i*3+channel]);
            }
        }

        if(DegreeSH >= 2){
            prevElems += elemsDegree*3;
            elemsDegree = 5;
            for(int i=0;i<elemsDegree;i++)
            {
                for(int channel=0;channel<3;channel++)
                    channelSH[channel].push_back(vecSH[prevElems+i*3+channel]);
            }
        }

        if(DegreeSH >= 3){
            prevElems += elemsDegree*3;
            uint curIdx = prevElems;

            while(curIdx < vecSH.size())
            {
                for(int channel=0;channel<3;channel++)
                    channelSH[channel].push_back(vecSH[curIdx+channel]);
                curIdx += 3;
            }
        }

        return channelSH;
    }

    static vector<ScalarType> mergeSHChannels(vector<vector<ScalarType>> &channelSH)
    {
        vector<ScalarType> vecSH;
        int elemsDegree, prevElems;

        // SH0
        vecSH.push_back(channelSH[0][0]);
        vecSH.push_back(channelSH[1][0]);
        vecSH.push_back(channelSH[2][0]);

        if(DegreeSH >= 1){
            prevElems = 1;
            elemsDegree = 3;
            for(int i=0;i<elemsDegree;i++)
            {
                for(int channel=0;channel<3;channel++)
                    vecSH.push_back(channelSH[channel][prevElems+i]);
            }
        }

        if(DegreeSH >= 2){
            prevElems += elemsDegree;
            elemsDegree = 5;
            for(int i=0;i<elemsDegree;i++)
            {
                for(int channel=0;channel<3;channel++)
                    vecSH.push_back(channelSH[channel][prevElems+i]);
            }
        }

        if(DegreeSH >= 3){
            prevElems += elemsDegree;
            uint curIdx = prevElems;

            while(curIdx < channelSH[0].size())
            {
                for(int channel=0;channel<3;channel++)
                    vecSH.push_back(channelSH[channel][curIdx]);
                curIdx += 1;
            }
        }

        return vecSH;
    }

public:
    constexpr static const float SH_C0 = 0.28209479177387814;

    // Constructors
    GaussianSplat() {}

    GaussianSplat(vcg::Quaternion<ScalarType> &rot, vcg::Point3<ScalarType> &scale,
                  vector<ScalarType> &vecSH, ScalarType opacity) {
        // Set rotation and scale
        this->rot = rot;
        this->scale = scale;

        vector<vector<ScalarType>> channelSH = splitSHChannels(vecSH);

        // Need to add 1 to degree to use Spherical Harmonics correctly
        this->sphR = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&channelSH[0][0]);  // red
        this->sphG = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&channelSH[1][0]);  // green
        this->sphB = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&channelSH[2][0]);  // blue

        // Save original values
        this->colorValues = vcg::Point4<ScalarType>(
            channelSH[0][0],
            channelSH[1][0],
            channelSH[2][0],
            opacity);

        // Set color
        this->color = vcg::Color4b(
            fdcToColor(SH_C0*channelSH[0][0]),
            fdcToColor(SH_C0*channelSH[1][0]),
            fdcToColor(SH_C0*channelSH[2][0]),
            opacityToAlpha(opacity)
        );
    }

    GaussianSplat(vcg::Quaternion<ScalarType> &rot, vcg::Point3<ScalarType> &scale, vcg::Color4b &color) {
        this->rot = rot;
        this->scale = scale;
        this->color = color;

        vector<ScalarType> shr = {colorToFdc(color[0])/SH_C0};
        vector<ScalarType> shg = {colorToFdc(color[1])/SH_C0};
        vector<ScalarType> shb = {colorToFdc(color[2])/SH_C0};

        // Need to add 1 to degree to use Spherical Harmonics correctly
        this->sphR = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&shr[0]);  // red
        this->sphG = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&shg[0]);  // green
        this->sphB = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&shb[0]);  // blue

        // Approximate original values
        this->colorValues = vcg::Point4<ScalarType>(
            shr[0],
            shg[0],
            shb[0],
            alphaToOpacity(color[3]));
    }

    GaussianSplat(vcg::Quaternion<ScalarType> &rot, vcg::Point3<ScalarType> &scale,
                  vcg::Color4b &color, vector<ScalarType> &vecSH) {
        // Set rotation and scale
        this->rot = rot;
        this->scale = scale;

        vector<vector<ScalarType>> channelSH = splitSHChannels(vecSH);
        channelSH[0].insert(channelSH[0].begin(), colorToFdc(color[0])/SH_C0);
        channelSH[1].insert(channelSH[1].begin(), colorToFdc(color[1])/SH_C0);
        channelSH[2].insert(channelSH[2].begin(), colorToFdc(color[2])/SH_C0);

        // Need to add 1 to degree to use Spherical Harmonics correctly
        this->sphR = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&channelSH[0][0]);  // red
        this->sphG = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&channelSH[1][0]);  // green
        this->sphB = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::Wrap(&channelSH[2][0]);  // blue

        // Save original values
        this->colorValues = vcg::Point4<ScalarType>(
            channelSH[0][0],
            channelSH[1][0],
            channelSH[2][0],
            color[3]);

        // Set color
        this->color = color;
    }

    ~GaussianSplat() {
    }

    static int getNumElemChannel()
    {
        int ret;

        if(DegreeSH == 0)
            ret = 1;
        else if(DegreeSH == 1)
            ret = 4;
        else if(DegreeSH == 2)
            ret = 9;
        else if(DegreeSH == 3)
            ret = 16;
        return ret;
    }

    static vector<ScalarType> getVectorSH(GaussianSplat &gs)
    {
        int numElems = getNumElemChannel();
        ScalarType* coeffR = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::getCoefficients(gs.sphR);
        ScalarType* coeffG = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::getCoefficients(gs.sphG);
        ScalarType* coeffB = vcg::math::SphericalHarmonics<ScalarType,DegreeSH+1>::getCoefficients(gs.sphB);

        vector<vector<ScalarType>> channelSH = {vector<ScalarType>(coeffR, coeffR+numElems), vector<ScalarType>(coeffG, coeffG+numElems), vector<ScalarType>(coeffB, coeffB+numElems)};

        return mergeSHChannels(channelSH);
    }

    static vcg::Point4<ScalarType> getColorValues(GaussianSplat &gs)
    {
        return gs.colorValues;
    }

    static vcg::Color4b getColor(GaussianSplat &gs)
    {
        return gs.color;
    }

    static vcg::Color4b getColor(GaussianSplat &gs, ScalarType theta, ScalarType phi)
    {
        return vcg::Color4b(fdcToColor(gs.sphR(theta, phi)), // r
                            fdcToColor(gs.sphG(theta, phi)), // g
                            fdcToColor(gs.sphB(theta, phi)), // b
                            gs.color[3]); // alpha
    }

    static vcg::Quaternion<ScalarType> getRotation(GaussianSplat &gs)
    {
        return gs.rot;
    }

    static vcg::Point3<ScalarType> getScale(GaussianSplat &gs)
    {
        return gs.scale;
    }
};

#endif // GAUSSIAN_SPLAT_H