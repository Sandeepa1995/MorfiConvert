#ifndef MORFICONVERT_FILE2OBJ_H
#define MORFICONVERT_FILE2OBJ_H

#include "GeneData.h"
#include <vector>
#include <iostream>

namespace File2Obj {
    //The NCBI BLAST data file conversion function. Supports .txt, .json and .xml formats.
    std::vector<GeneData> getInNCBIBLASTData(std::string path,int thrdCnt=10);

    // Get custom BLAST data in a custom file (not NCBI)
    std::vector<GeneData> getInCustomBLASTData(std::string path, int thrdCnt=10);

    // Get custom MSA(Multiple Sequence Alignment) data in a custom file (not NCBI)
    std::vector<GeneData> getInCustomMSAData(std::string path);

    //The FASTA data file conversion function
    std::vector<GeneData> getInFASTAData(std::string path,int thrdCnt=10);

    //The GCG data file conversion function.
    std::vector<GeneData> getInGCGData(std::string path,int thrdCnt=10);

    //The EMBL data file conversion function.
    std::vector<GeneData> getInEMBLData(std::string path,int thrdCnt=10);

    //The Genbank data file conversion function.
    std::vector<GeneData> getInGenebankData(std::string path,int thrdCnt=10);

    //The PIR data file conversion function.
    std::vector<GeneData> getInPIRData(std::string path,int thrdCnt=10);

    //The PHYLI data file conversion function.
    std::vector<GeneData> getInPHYLIPData(std::string path,int thrdCnt=10);

    //The Clustral data file (A Multiple Sequence Alignment) conversion function.
    std::vector<GeneData> getInClustralData(std::string path);

    //The MSF data file conversion function.
    std::vector<GeneData> getInMSFData(std::string path);

    //The FASTA Report data file conversion function.
    std::vector<GeneData> getInFASTARepData(std::string path,int thrdCnt=10);
}

#endif