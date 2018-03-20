#ifndef MORFICONVERT_MORFICONVERT_H
#define MORFICONVERT_MORFICONVERT_H

#include <iostream>
#include <vector>

#include "GeneData.h"

namespace morfiConvert{
    //Components / filters
    std::vector<GeneData> getInData(std::string path, std::string origin);
    std::vector<std::string> giveOutData(std::vector<GeneData> res, std::string outType);
    void writeOutData(std::vector<GeneData> res, std::string path, std::string outType);

    //Intermediary
    std::string identify(std::string path);
    void configLocal(std::string path);
    std::string completeGene(std::string ncbiID, std::string path);
    std::vector<GeneData> completeGenes(std::vector<GeneData> geneSet);
    std::vector<GeneData> multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex);
    std::vector<GeneData> multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString);
    GeneData indexFragment(GeneData inpData, int strtIndex, int endIndex);
    GeneData stringFragment(GeneData inpData, std::string strtString, std::string endString);

//    void MSAConvert(std::string inPath);

    //Pipeline
    std::vector<std::string> objConvert(std::string inPath, std::string origin, std::string outType);
    void fileConvert(std::string inPath, std::string outPath, std::string origin, std::string outType);
}

#endif