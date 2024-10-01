#ifndef IMPORT_SPLAT_H
#define IMPORT_SPLAT_H

#include <iostream>
#include <fstream>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>
#include "./gaussian_splat.h"

namespace vcg {
namespace tri {
namespace io {

template <class OpenMeshType, int DegreeSH>
class ImporterSPLAT
{
public:
    typedef typename OpenMeshType::VertexType VertexType;
    typedef typename OpenMeshType::VertexIterator GaussianIterator;

    static int Open(OpenMeshType &m,  const char * filename, CallBackPos *cb=0)
    {
        assert(filename != 0);
        assert(DegreeSH >= 0 && DegreeSH <= 3);

        const int byteRowSize = 32;

        int j;
        GaussianIterator gi;

        auto handleGS = vcg::tri::Allocator<OpenMeshType>:: template AddPerVertexAttribute<GaussianSplat<float,DegreeSH>> (m, "gs");

        if(!tri::Allocator<OpenMeshType>::IsValidHandle(m, handleGS))
        {
            return vcg::ply::E_ELEMNOTFOUND;
        }

        ifstream inputFile(filename, std::ios::in | std::ios::binary);
        if (!inputFile)
        {
            return vcg::ply::E_CANTOPEN;
        }

        // copies all data into buffer
        std::vector<unsigned char> buffer(std::istreambuf_iterator<char>(inputFile), {});

        // Add vertices
        int splatCount = buffer.size()/byteRowSize;
        vcg::tri::Allocator<OpenMeshType>::AddVertices(m, splatCount);

        // Read a total of 32 bytes from buffer
        Point3f point; // 3* 4 bytes
        Point3f scale; // 3* 4 bytes
        Point4<unsigned char> tmpColor; // 4 bytes
        Point4<unsigned char> tmpRot; // 4 bytes

        for(j=0,gi=m.vert.begin();gi!=m.vert.end();++gi,j++){
            //((m.vn+m.fn) != 0) all vertices and faces have been marked as deleted but the are still in the vert/face vectors
            if(cb && ((j%1000)==0) && ((m.vn) != 0) )(*cb)( (100*j)/(m.vn), "Loading Gaussians");

            memcpy(&point[0], &buffer[0] + j*byteRowSize, sizeof(float)*3);
            memcpy(&scale[0], &buffer[0] + j*byteRowSize + sizeof(float)*3, sizeof(float)*3);
            memcpy(&tmpColor[0], &buffer[0] + j*byteRowSize + sizeof(float)*6, sizeof(unsigned char)*4);
            memcpy(&tmpRot[0], &buffer[0] + j*byteRowSize + sizeof(float)*6 + sizeof(unsigned char)*4, sizeof(unsigned char)*4);

            gi->P() = point;
            Color4b color(tmpColor);
            Quaternion<float> rot((float(tmpRot[0]) - 128) / 128, (float(tmpRot[1]) - 128) / 128, (float(tmpRot[2]) - 128) / 128, (float(tmpRot[3]) - 128) / 128);

            handleGS[gi] = GaussianSplat<float,DegreeSH>(rot, scale, color);
        }

        return vcg::ply::E_NOERROR;
    }
};

}
}
}

#endif // IMPORT_SPLAT_H
