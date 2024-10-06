#ifndef IMPORT_KSPLAT_H
#define IMPORT_KSPLAT_H

#include <iostream>
#include <fstream>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>
#include "./gaussian_splat.h"

#include <Eigen/src/Core/arch/CUDA/Half.h>

namespace vcg {
namespace tri {
namespace io {

template <class OpenMeshType, int DegreeSH>
class ImporterKSPLAT
{
private:
    // KSPLAT organisation: main header -> all section headers -> all section data, each with all the buckets at the start
    static const int mainHeaderSizeBytes = 4096;
    static const int sectionHeaderSizeBytes = 1024;
    constexpr static float defaultSphericalHarmonics8BitCompressionRange = 3;
    constexpr static int sphericalHarmonicsComponentCountForDegree[] = {0, 9, 24}; // Limited to SH < 3

    struct compressionLevelInfo
    {
        int bytesPerCenter;
        int bytesPerScale;
        int bytesPerRotation;
        int bytesPerColor;
        int scaleOffsetBytes;
        int rotationOffsetBytes;
        int colorOffsetBytes;
        int sphericalHarmonicsOffsetBytes;
        int scaleRange;
        int bytesPerSphericalHarmonicsComponent;
        int sphericalHarmonicsOffsetFloat;
        int sphericalHarmonicsDegrees[3];
    };

    constexpr static compressionLevelInfo compressionInfo[3] = {
        {12, 12, 16, 4, 12, 24, 40, 44, 1, 4, 11, {44, 80, 140}},
        {6, 6, 8, 4, 6, 12, 20, 24, 32767, 2, 12, {24, 42, 72}},
        {6, 6, 8, 4, 6, 12, 20, 24, 32767, 1, 12, {24, 33, 48}}
    };

    struct mainHeader
    {
        uint8_t majorVersion;
        uint8_t minorVersion;
        uint32_t maxSectionCount;
        uint32_t sectionCount;
        uint32_t maxSplatCount;
        uint32_t splatCount;
        uint16_t compressionLevel;
        float sceneCenter[3];
        float minSphericalHarmonicsCoeff;
        float maxSphericalHarmonicsCoeff;
    };

    static int ReadMainHeader(uint8_t *buffer, mainHeader &header)
    {
        int curIdx = 0;

        memcpy(&header.majorVersion, &buffer[0], sizeof(uint8_t));
        curIdx += sizeof(uint8_t);
        memcpy(&header.minorVersion, &buffer[0] + curIdx, sizeof(uint8_t));
        curIdx += sizeof(uint8_t) * 3; // Need to account for unused bytes
        memcpy(&header.maxSectionCount, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.sectionCount, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.maxSplatCount, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.splatCount, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.compressionLevel, &buffer[0] + curIdx, sizeof(uint16_t));
        curIdx += sizeof(uint16_t) * 2; // alignment
        memcpy(&(header.sceneCenter[0]), &buffer[0] + curIdx, sizeof(float)*3);
        curIdx += sizeof(float) * 3;
        memcpy(&header.minSphericalHarmonicsCoeff, &buffer[0] + curIdx, sizeof(float));
        curIdx += sizeof(float);
        memcpy(&header.maxSphericalHarmonicsCoeff, &buffer[0] + curIdx, sizeof(float));

        return 0;
    }

    struct sectionHeader
    {
        // Main vars (written and read from file)
        uint32_t splatCount;
        uint32_t maxSplatCount;
        uint32_t bucketSize;
        uint32_t bucketCount;
        float bucketBlockSize;
        uint16_t bucketStorageSizeBytes;
        uint32_t compressionScaleRange;
        uint32_t storageSizeBytes;
        uint32_t fullBucketCount;
        uint32_t partiallyFilledBucketCount;
        uint16_t sphericalHarmonicsDegree;
/*
        // Derived vars
        bytesPerSplat;
        splatCountOffset;
        halfBucketBlockSize;
        bucketsStorageSizeBytes;
        splatDataStorageSizeBytes;
        compressionScaleFactor;
        base;
        bucketsBase;
        dataBase;*/
    };

    static int ReadSectionHeader(uint8_t *buffer, sectionHeader &header)
    {
        int curIdx = 0;

        memcpy(&header.splatCount, &buffer[0], sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.maxSplatCount, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.bucketSize, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.bucketCount, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.bucketBlockSize, &buffer[0] + curIdx, sizeof(float));
        curIdx += sizeof(float);
        memcpy(&header.bucketStorageSizeBytes, &buffer[0] + curIdx, sizeof(uint16_t));
        curIdx += sizeof(uint16_t) * 2; // alignment
        memcpy(&header.compressionScaleRange, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.storageSizeBytes, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.fullBucketCount, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.partiallyFilledBucketCount, &buffer[0] + curIdx, sizeof(uint32_t));
        curIdx += sizeof(uint32_t);
        memcpy(&header.sphericalHarmonicsDegree, &buffer[0] + curIdx, sizeof(uint16_t));

        return 0;
    }


public:
    typedef typename OpenMeshType::VertexType VertexType;
    typedef typename OpenMeshType::VertexIterator GaussianIterator;

    static int Open(OpenMeshType &m,  const char * filename, CallBackPos *cb=0)
    {
        assert(filename != 0);
        assert(DegreeSH >= 0 && DegreeSH <= 3);

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

        mainHeader header;
        ReadMainHeader(&buffer[0], header);

        if (header.majorVersion > 0 || (header.majorVersion == 0 && header.minorVersion > 1))
        {
            return vcg::ply::E_CANTOPEN;
        }

        // Add vertices
        vcg::tri::Allocator<OpenMeshType>::AddVertices(m, header.splatCount);

        int bytesPerCenter = compressionInfo[header.compressionLevel].bytesPerCenter;
        int bytesPerScale = compressionInfo[header.compressionLevel].bytesPerScale;
        int bytesPerRotation = compressionInfo[header.compressionLevel].bytesPerRotation;
        int bytesPerColor = compressionInfo[header.compressionLevel].bytesPerColor;

        int sectionHeaderBase = mainHeaderSizeBytes;
        int sectionBase = sectionHeaderBase + header.maxSectionCount * sectionHeaderSizeBytes;
        int bufferOffset = 0;
        int vertIdx = 0;

        for(uint32_t i=0;i<header.sectionCount;i++){
            sectionHeader curHeader;
            ReadSectionHeader(&buffer[0] + sectionHeaderBase + sectionHeaderSizeBytes*i, curHeader);

            int compressionScaleRange = curHeader.compressionScaleRange;
            if(compressionScaleRange == 0)
                compressionScaleRange = compressionInfo[header.compressionLevel].scaleRange;

            float compressionScaleFactor = (curHeader.bucketBlockSize/2) / compressionScaleRange;

            int sphericalHarmonicsComponentsPerSplat = sphericalHarmonicsComponentCountForDegree[curHeader.sphericalHarmonicsDegree];
            int sphericalHarmonicsBytesPerSplat = compressionInfo[header.compressionLevel].bytesPerSphericalHarmonicsComponent * sphericalHarmonicsComponentsPerSplat;

            int bytesPerSplat = bytesPerCenter +
                                bytesPerScale +
                                bytesPerRotation +
                                bytesPerColor +
                                sphericalHarmonicsBytesPerSplat;

            // Bucket-related vars
            int bucketsMetaDataSizeBytes = header.compressionLevel >= 1 ? curHeader.partiallyFilledBucketCount * sizeof(uint32_t) : 0;
            int bucketData = header.compressionLevel >= 1 ? bucketsMetaDataSizeBytes + curHeader.bucketStorageSizeBytes * curHeader.bucketCount : 0;

            for(uint32_t j=0;j<curHeader.splatCount;j++){

                if(cb && ((vertIdx%1000)==0) && ((m.vn) != 0) )(*cb)( (100*vertIdx)/(m.vn), "Loading Gaussians");
                bufferOffset = j*bytesPerSplat;

                // Position
                Point3f pos;

                if(header.compressionLevel == 0){
                    memcpy(&pos[0], &buffer[0] + sectionBase + bufferOffset, sizeof(float) * 3);
                } else {
                    uint32_t maxSplatIndexInFullBuckets = curHeader.fullBucketCount * curHeader.bucketSize;
                    int bucketIdx = floor(j/curHeader.bucketSize);

                    if(j >= maxSplatIndexInFullBuckets) {
                        uint32_t bucketSplatIndex = maxSplatIndexInFullBuckets;
                        bucketIdx = curHeader.fullBucketCount;
                        uint32_t partiallyFullBucketIndex = 0;
                        vector<uint32_t> partiallyFilledBucketLengths(curHeader.partiallyFilledBucketCount);

                        if(curHeader.partiallyFilledBucketCount>0)
                            memcpy(&partiallyFilledBucketLengths[0], &buffer[0] + sectionBase, sizeof(uint32_t) * curHeader.partiallyFilledBucketCount);

                        while (bucketSplatIndex < curHeader.splatCount) {
                            uint32_t currentPartiallyFilledBucketSize = partiallyFilledBucketLengths[partiallyFullBucketIndex];
                            if (j >= bucketSplatIndex && j < bucketSplatIndex + currentPartiallyFilledBucketSize) {
                                break;
                            }
                            bucketSplatIndex += currentPartiallyFilledBucketSize;
                            bucketIdx++;
                            partiallyFullBucketIndex++;
                        }
                    }

                    int bucketStorageSizeFloats = 3;
                    int bucketOffset = bucketsMetaDataSizeBytes + bucketIdx * bucketStorageSizeFloats * sizeof(float);
                    Point3f bucket;
                    uint16_t tmpPos[3]; // Need to read compressed float16

                    memcpy(&tmpPos[0], &buffer[0] + sectionBase + bucketData + bufferOffset, sizeof(uint16_t) * 3);
                    pos = Point3f(tmpPos[0], tmpPos[1], tmpPos[2]);

                    memcpy(&bucket[0], &buffer[0] + sectionBase + bucketOffset, sizeof(float) * bucketStorageSizeFloats);
                    pos = (pos - Point3f(compressionScaleRange, compressionScaleRange, compressionScaleRange)) * compressionScaleFactor + bucket;
                }

                // Scale
                Point3f scale;
                if(header.compressionLevel > 0){
                    Eigen::half tmpScale[3];
                    memcpy(&tmpScale[0], &buffer[0] + sectionBase + bucketData + bufferOffset + compressionInfo[header.compressionLevel].scaleOffsetBytes, sizeof(uint16_t) * 3);
                    scale = Point3f(Eigen::half_impl::half_to_float(tmpScale[0]), Eigen::half_impl::half_to_float(tmpScale[1]), Eigen::half_impl::half_to_float(tmpScale[2]));
                } else {
                    memcpy(&scale[0], &buffer[0] + sectionBase + bucketData + bufferOffset + compressionInfo[header.compressionLevel].scaleOffsetBytes, sizeof(float) * 3);
                }

                // Rotation
                Quaternion<float> rot;
                if(header.compressionLevel > 0){
                    Eigen::half tmpRot[4];
                    memcpy(&tmpRot[0], &buffer[0] + sectionBase + bucketData + bufferOffset + compressionInfo[header.compressionLevel].rotationOffsetBytes, sizeof(uint16_t) * 4);
                    rot = Quaternion<float>(Eigen::half_impl::half_to_float(tmpRot[0]), Eigen::half_impl::half_to_float(tmpRot[1]),
                                            Eigen::half_impl::half_to_float(tmpRot[2]), Eigen::half_impl::half_to_float(tmpRot[3]));
                } else {
                    memcpy(&rot[0], &buffer[0] + sectionBase + bucketData + bufferOffset + compressionInfo[header.compressionLevel].rotationOffsetBytes, sizeof(float) * 4);
                }

                // Color
                Color4b color;
                memcpy(&color[0], &buffer[0] + sectionBase + bucketData + bufferOffset + compressionInfo[header.compressionLevel].colorOffsetBytes, sizeof(uint8_t) * 4);

                // Spherical Harmonics
                if(curHeader.sphericalHarmonicsDegree > 0){
                    int numSh = sphericalHarmonicsComponentCountForDegree[curHeader.sphericalHarmonicsDegree];
                    vector<float> shVec(numSh);
                    int shStart = sectionBase + bucketData + bufferOffset + compressionInfo[header.compressionLevel].sphericalHarmonicsOffsetBytes;

                    if(header.compressionLevel == 0){
                        memcpy(&shVec[0], &buffer[0] + shStart, sizeof(float) * numSh);
                    } else if(header.compressionLevel == 1){
                        vector<Eigen::half> tmpSh(numSh);
                        memcpy(&tmpSh[0], &buffer[0] + shStart, sizeof(uint16_t) * numSh);
                        for(int k=0;k<numSh;k++){
                            shVec[k] = Eigen::half_impl::half_to_float(tmpSh[k]);
                        }
                    } else if(header.compressionLevel == 2){
                        // TODO: implement decompression
                        vector<uint8_t> tmpSh8bit(numSh);
                        memcpy(&shVec[0], &buffer[0] + shStart, sizeof(uint8_t) * numSh);
                        for(int k=0;k<numSh;k++){
                            float range = header.maxSphericalHarmonicsCoeff - header.minSphericalHarmonicsCoeff;
                            shVec[k] = tmpSh8bit[k] / 255 * range + header.minSphericalHarmonicsCoeff;
                        }
                    }

                    handleGS[vertIdx] = GaussianSplat<float,DegreeSH>(rot, scale, color, shVec);
                } else {
                    handleGS[vertIdx] = GaussianSplat<float,DegreeSH>(rot, scale, color);
                }

                m.vert[vertIdx].P() = pos;

                vertIdx++;
            }
            sectionBase += bucketData + bufferOffset;
        }

        return vcg::ply::E_NOERROR;
    }
};

}
}
}

#endif // IMPORT_KSPLAT_H
