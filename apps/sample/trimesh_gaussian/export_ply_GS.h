#ifndef EXPORT_PLY_GS_H
#define EXPORT_PLY_GS_H

#include<stddef.h>
#include<cmath>
#include<wrap/callback.h>
#include<wrap/ply/plylib.h>
#include<wrap/io_trimesh/io_mask.h>
#include<wrap/io_trimesh/io_ply.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<wrap/io_trimesh/export_ply.h>
#include"./gaussian_splat.h"


namespace vcg {
namespace tri {
namespace io {


template <class SaveMeshType, int DegreeSH>
class ExporterPLYGS
{
private:
    static int clamp(int v, int lo, int hi) {
        return min(hi, max(lo, v));
    }
public:
    static int findPropertyIdx(const vector<string> properties, int propLen, const std::string name)
    {
        return find(&properties[0], &properties[0] + propLen, name) - &properties[0];
    }

    static int Save(SaveMeshType &m, const char * filename, bool binary=true)
    {
        PlyInfo pi;
        return Save(m,filename,binary,pi);
    }

    static int Save(SaveMeshType &m,  const char * filename, int savemask, bool binary = true, CallBackPos *cb=0)
    {
        PlyInfo pi;
        pi.mask=savemask;
        return Save(m,filename,binary,pi,cb);
    }

    // Save a mesh with all the valid ply attributes for Gaussian Splat objects, returns 0 on success.
    static int Save(SaveMeshType &m,  const char * filename, bool binary, PlyInfo &pi, CallBackPos *cb=0)
    {
        assert(filename != 0);
        assert(DegreeSH >= 0 && DegreeSH <= 3);

        if(!vcg::tri::HasPerVertexAttribute(m,"gs"))
        {
            return ply::E_PROPNOTFOUND;
        }

        const vector<string> properties = {"f_dc_0", "f_dc_1", "f_dc_2", "f_rest_0", "f_rest_1", "f_rest_2", "f_rest_3", "f_rest_4", "f_rest_5", "f_rest_6", "f_rest_7", "f_rest_8", "f_rest_9", "f_rest_10", "f_rest_11", "f_rest_12", "f_rest_13", "f_rest_14", "f_rest_15", "f_rest_16", "f_rest_17", "f_rest_18", "f_rest_19", "f_rest_20", "f_rest_21", "f_rest_22", "f_rest_23", "f_rest_24", "f_rest_25", "f_rest_26", "f_rest_27", "f_rest_28", "f_rest_29", "f_rest_30", "f_rest_31", "f_rest_32", "f_rest_33", "f_rest_34", "f_rest_35", "f_rest_36", "f_rest_37", "f_rest_38", "f_rest_39", "f_rest_40", "f_rest_41", "f_rest_42", "f_rest_43", "f_rest_44", "opacity", "scale_0", "scale_1", "scale_2", "rot_0", "rot_1", "rot_2", "rot_3"};
        typename SaveMeshType::template PerVertexAttributeHandle<GaussianSplat<float,DegreeSH>> handleGS =
            vcg::tri::Allocator<SaveMeshType>:: template GetPerVertexAttribute<GaussianSplat<float,DegreeSH>> (m, "gs");
        std::vector<typename SaveMeshType::template PerVertexAttributeHandle<float>> handleVec(properties.size());

        // Save normals using pi mask
        if( vcg::tri::HasPerVertexNormal(m) )
            pi.mask |= Mask::IOM_VERTNORMAL;

        // Init PlyInfo with custom attributes
        for(uint i=0;i<properties.size();i++) {

            // Add attributes to PlyInfo
            pi.AddPerVertexFloatAttribute(properties[i]);

            // Add vertex attributes if missing
            if (!vcg::tri::HasPerVertexAttribute(m, properties[i])) {
                handleVec[i] = vcg::tri::Allocator<SaveMeshType>::template AddPerVertexAttribute<float>(m, properties[i]);
            } else {
                handleVec[i] = vcg::tri::Allocator<SaveMeshType>::template GetPerVertexAttribute<float>(m, properties[i]);
            }
        }

        // Find indices of relevant properties
        // Color (technically sh0)
        int propFDC0 = findPropertyIdx(properties, properties.size(), "f_dc_0");
        int propFDC1 = findPropertyIdx(properties, properties.size(), "f_dc_1");
        int propFDC2 = findPropertyIdx(properties, properties.size(), "f_dc_2");
        int propOpacity = findPropertyIdx(properties, properties.size(), "opacity");
        // Scale
        int propScaleX = findPropertyIdx(properties, properties.size(), "scale_0");
        int propScaleY = findPropertyIdx(properties, properties.size(), "scale_1");
        int propScaleZ = findPropertyIdx(properties, properties.size(), "scale_2");
        // Rotation
        int propRot0 = findPropertyIdx(properties, properties.size(), "rot_0");
        int propRot1 = findPropertyIdx(properties, properties.size(), "rot_1");
        int propRot2 = findPropertyIdx(properties, properties.size(), "rot_2");
        int propRot3 = findPropertyIdx(properties, properties.size(), "rot_3");
        // Spherical Harmonics: only first and last, assume ordered list
        int propStartSH = findPropertyIdx(properties, properties.size(), "f_rest_0");
        int propEndSH = findPropertyIdx(properties, properties.size(), "f_rest_44");

        // Loop through vertices, where each one represents a gaussian
        // Reference: https://github.com/limacv/GaussianSplattingViewer/blob/main/shaders/gau_vert.glsl
        for(typename SaveMeshType::VertexIterator gi=m.vert.begin();gi!=m.vert.end();++gi)
        {
            if(cb && (vcg::tri::Index(m, *gi)%1000==0) && (m.vn != 0) )(*cb)(100*tri::Index(m, *gi)/m.vn, "Saving Splats");

            GaussianSplat<float,DegreeSH> gs = handleGS[gi];

            // Set scale
            handleVec[propScaleX][gi] = log(GaussianSplat<float,DegreeSH>::getScale(gs)[0]);
            handleVec[propScaleY][gi] = log(GaussianSplat<float,DegreeSH>::getScale(gs)[1]);
            handleVec[propScaleZ][gi] = log(GaussianSplat<float,DegreeSH>::getScale(gs)[2]);

            // Set rotation quaternion
            handleVec[propRot0][gi] = GaussianSplat<float,DegreeSH>::getRotation(gs)[0];
            handleVec[propRot1][gi] = GaussianSplat<float,DegreeSH>::getRotation(gs)[1];
            handleVec[propRot2][gi] = GaussianSplat<float,DegreeSH>::getRotation(gs)[2];
            handleVec[propRot3][gi] = GaussianSplat<float,DegreeSH>::getRotation(gs)[3];

            // Set Spherical Harmonics coefficients
            std::vector<float> vecSH = GaussianSplat<float,DegreeSH>::getVectorSH(gs);

            handleVec[propFDC0][gi] = vecSH[0];
            handleVec[propFDC1][gi] = vecSH[1];
            handleVec[propFDC2][gi] = vecSH[2];

            for(int i=propStartSH;i<=propEndSH;i++)
            {
                handleVec[i][gi] = vecSH[3+(i-propStartSH)];
            }

            // Set alpha value
            handleVec[propOpacity][gi] = GaussianSplat<float,DegreeSH>::getColorValues(gs)[3];

        }

        int ret = tri::io::ExporterPLY<SaveMeshType>::Save(m, filename, binary, pi);
        if(ret!=0)
        {
            printf("Unable to save %s for '%s'\n", filename, vcg::tri::io::ExporterPLY<SaveMeshType>::ErrorMsg(ret));
            return pi.status;
        }

        return 0;
    }

}; // end class


}
}
}
#endif // EXPORT_PLY_GS_H
