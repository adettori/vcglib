#ifndef EXPORT_SPLAT_H
#define EXPORT_SPLAT_H

#include <iostream>
#include <fstream>
#include <wrap/io_trimesh/export_ply.h>
#include "./gaussian_splat.h"

namespace vcg {
namespace tri {
namespace io {

template <class SaveMeshType, int DegreeSH>
class ExporterSPLAT
{
public:
    typedef typename SaveMeshType::VertexType VertexType;
    typedef typename SaveMeshType::VertexIterator GaussianIterator;

    static int Save(SaveMeshType &m,  const char * filename, CallBackPos *cb=0)
    {
        assert(filename != 0);
        assert(DegreeSH >= 0 && DegreeSH <= 3);

        const int byteRowSize = 32;
        int splatCount = m.vert.size();

        int j;
        GaussianIterator gi;

        auto handleGS = vcg::tri::Allocator<SaveMeshType>:: template GetPerVertexAttribute<GaussianSplat<float,DegreeSH>> (m, "gs");

        if(!tri::Allocator<SaveMeshType>::IsValidHandle(m, handleGS))
        {
            return vcg::ply::E_ELEMNOTFOUND;
        }

        // Init buffer that will be written to file
        std::vector<unsigned char> buffer(splatCount*byteRowSize, {});

        for(j=0,gi=m.vert.begin();gi!=m.vert.end();++gi,j++){
            //((m.vn+m.fn) != 0) all vertices and faces have been marked as deleted but the are still in the vert/face vectors
            if(cb && ((j%1000)==0) && ((m.vn) != 0) )(*cb)( (100*j)/(m.vn), "Saving Gaussians");

            Point3f scale = GaussianSplat<float,DegreeSH>::getScale(handleGS[gi]);
            Color4b color = GaussianSplat<float,DegreeSH>::getColor(handleGS[gi]);
            Quaternion<float> rotF = GaussianSplat<float,DegreeSH>::getRotation(handleGS[gi]);
            Quaternion<unsigned char> rotB((unsigned char) (rotF[0]*128+128), (unsigned char) (rotF[1]*128+128), (unsigned char) (rotF[2]*128+128), (unsigned char) (rotF[3]*128+128));

            memcpy(&buffer[0] + j*byteRowSize, &(gi->P()[0]), sizeof(float)*3);
            memcpy(&buffer[0] + j*byteRowSize + sizeof(float)*3, &(scale[0]), sizeof(float)*3);
            memcpy(&buffer[0] + j*byteRowSize + sizeof(float)*6, &(color[0]), sizeof(unsigned char)*4);
            memcpy(&buffer[0] + j*byteRowSize + sizeof(float)*6 + sizeof(unsigned char)*4, &(rotB[0]), sizeof(unsigned char)*4);
        }

        ofstream outputFile(filename, std::ios::out | std::ios::binary);
        if (!outputFile)
        {
            return vcg::ply::E_CANTOPEN;
        }
        std::ostream_iterator<unsigned char> outputIterator(outputFile);
        std::copy(buffer.begin(), buffer.end(), outputIterator);

        return vcg::ply::E_NOERROR;
    }
};

}
}
}

#endif // EXPORT_SPLAT_H
