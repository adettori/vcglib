/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2024                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef IMPORT_PLY_GS_H
#define IMPORT_PLY_GS_H

#include<stddef.h>
#include<wrap/callback.h>
#include<wrap/ply/plylib.h>
#include<wrap/io_trimesh/io_mask.h>
#include<wrap/io_trimesh/io_ply.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<wrap/io_trimesh/import_ply.h>
#include"./gaussian_splat.h"


namespace vcg {
namespace tri {
namespace io {

/**
This class encapsulate a filter for opening ply Gaussian Splat meshes.
The ply file format is quite extensible...
*/
template <class OpenMeshType, int DegreeSH>
class ImporterPLYGS
{
public:
    static int findPropertyIdx(const std::string properties[], int propLen, const std::string name)
    {
        return find(&properties[0], properties + propLen, name) - properties;
    }

    /// Standard call for reading a mesh, returns 0 on success.
    typename OpenMeshType::template PerVertexAttributeHandle<GaussianSplat<float,DegreeSH>>
    static Open( OpenMeshType &m, const char * filename, CallBackPos *cb=0)
    {
        PlyInfo pi;
        pi.cb=cb;
        return Open(m, filename, pi);
    }

    /// Read a mesh and store in loadmask the loaded field
    /// Note that loadmask is not read! just modified. You cannot specify what fields
    /// have to be read. ALL the data for which your mesh HasSomething and are present
    /// in the file are read in.
    static int Open( OpenMeshType &m, const char * filename, int & loadmask, CallBackPos *cb =0)
    {
        PlyInfo pi;
        pi.cb=cb;
        int r = Open(m, filename,pi);
        loadmask=pi.mask;
        return r;
    }

    /// read a mesh with all the possible option specified in the PlyInfo obj, returns 0 on success.
    static int Open( OpenMeshType &m, const char * filename, PlyInfo &pi)
    {
        assert(filename != 0);
        assert(DegreeSH >= 0 && DegreeSH <= 3);
        const std::string properties[] = {"f_dc_0", "f_dc_1", "f_dc_2", "f_rest_0", "f_rest_1", "f_rest_2", "f_rest_3", "f_rest_4", "f_rest_5", "f_rest_6", "f_rest_7", "f_rest_8", "f_rest_9", "f_rest_10", "f_rest_11", "f_rest_12", "f_rest_13", "f_rest_14", "f_rest_15", "f_rest_16", "f_rest_17", "f_rest_18", "f_rest_19", "f_rest_20", "f_rest_21", "f_rest_22", "f_rest_23", "f_rest_24", "f_rest_25", "f_rest_26", "f_rest_27", "f_rest_28", "f_rest_29", "f_rest_30", "f_rest_31", "f_rest_32", "f_rest_33", "f_rest_34", "f_rest_35", "f_rest_36", "f_rest_37", "f_rest_38", "f_rest_39", "f_rest_40", "f_rest_41", "f_rest_42", "f_rest_43", "f_rest_44", "opacity", "scale_0", "scale_1", "scale_2", "rot_0", "rot_1", "rot_2", "rot_3"};

        const int propLen = sizeof(properties)/sizeof(properties[0]);
        std::vector<typename OpenMeshType::template PerVertexAttributeHandle<float>> handleVec(propLen);
        auto handleGs = vcg::tri::Allocator<OpenMeshType>:: template AddPerVertexAttribute<GaussianSplat<float,DegreeSH>> (m, "gs");

        // Init attribute handler
        for(int i=0;i<propLen;i++) {
            pi.AddPerVertexFloatAttribute(properties[i]);
            // add a per-vertex attribute with type float
            handleVec[i] = vcg::tri::Allocator<OpenMeshType>:: template GetPerVertexAttribute<float> (m, properties[i]);
        }

        // Load gs ply file
        int ret = vcg::tri::io::ImporterPLY<OpenMeshType>::Open(m, filename, pi);
        if(ret!=0)
        {
            printf("Unable to open %s for '%s'\n", filename, vcg::tri::io::ImporterPLY<OpenMeshType>::ErrorMsg(ret));
            return pi.status;
        }

        // Find indices of relevant properties
        // Color (technically sh0)
        int propFDC0 = findPropertyIdx(properties, propLen, "f_dc_0");
        int propFDC1 = findPropertyIdx(properties, propLen, "f_dc_1");
        int propFDC2 = findPropertyIdx(properties, propLen, "f_dc_2");
        int propOpacity = findPropertyIdx(properties, propLen, "opacity");
        // Scale
        int propScaleX = findPropertyIdx(properties, propLen, "scale_0");
        int propScaleY = findPropertyIdx(properties, propLen, "scale_1");
        int propScaleZ = findPropertyIdx(properties, propLen, "scale_2");
        // Rotation
        int propRot0 = findPropertyIdx(properties, propLen, "rot_0");
        int propRot1 = findPropertyIdx(properties, propLen, "rot_1");
        int propRot2 = findPropertyIdx(properties, propLen, "rot_2");
        int propRot3 = findPropertyIdx(properties, propLen, "rot_3");
        // Spherical Harmonics: only first, assume ordered list
        int propStartSH = findPropertyIdx(properties, propLen, "f_rest_0");
        int propEndSH = findPropertyIdx(properties, propLen, "f_rest_44");

        // Loop through vertices, where each one represents a gaussian
        // Reference: https://github.com/limacv/GaussianSplattingViewer/blob/main/shaders/gau_vert.glsl
        for(typename OpenMeshType::VertexIterator gi=m.vert.begin();gi!=m.vert.end();++gi)
        {
            if(pi.cb && (vcg::tri::Index(m, *gi)%1000==0) && (m.vn != 0) )(*pi.cb)(100*tri::Index(m, *gi)/m.vn, "Loading Splats");

            vector<float> vecSH;
            // SH0
            vecSH.push_back(handleVec[propFDC0][gi]);
            vecSH.push_back(handleVec[propFDC1][gi]);
            vecSH.push_back(handleVec[propFDC2][gi]);

            // Everything else SH related
            for(int i=propStartSH;i<=propEndSH;i++)
            {
                vecSH.push_back(handleVec[i][gi]);
            }

            // Set rotation and scale
            // Need to exponentiate the scale values
            vcg::Quaternion<float> rotQuat = vcg::Quaternion<float>(handleVec[propRot0][gi], handleVec[propRot1][gi], handleVec[propRot2][gi], handleVec[propRot3][gi]);
            vcg::Point3<float> scale = vcg::Point3<float>(exp(handleVec[propScaleX][gi]), exp(handleVec[propScaleY][gi]), exp(handleVec[propScaleZ][gi]));
            handleGs[gi] = GaussianSplat<float,DegreeSH>(rotQuat, scale, vecSH, handleVec[propOpacity][gi]);
        }

        return 0;
    }

}; // end class


}
}
}
#endif // IMPORT_PLY_GS_H
