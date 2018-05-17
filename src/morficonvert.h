#ifndef MORFICONVERT_MORFICONVERT_H
#define MORFICONVERT_MORFICONVERT_H

#include <iostream>
#include <vector>

#include "GeneData.h"

namespace morfiConvert{
    //******************************************************************************************************************
    // Input components / filters / conversions
    //******************************************************************************************************************

    // Conversion of the Input Data to GeneData objects (text / character files)
    std::vector<GeneData> getInData(std::string path, std::string origin, int thrdCnt=10);

    // Conversion of the Input MSA Data to GeneData objects (Image files)
    std::vector<GeneData> getInImageData(std::vector<std::string> genes, std::string filePath);


    //******************************************************************************************************************
    // Output components / filters / conversions
    //******************************************************************************************************************

    // Conversion of GeneData objects to the requested string format
    std::vector<std::string> giveOutData(std::vector<GeneData> res, std::string outType, int thrdCnt=10);

    //Write the given GeneData objects to a specific file format
    void writeOutData(std::vector<GeneData> res, std::string path, std::string outType, int thrdCnt=10);


    //******************************************************************************************************************
    //Intermediary
    //******************************************************************************************************************

    //Identify type
    std::string identify(std::string path);

    //Configure files for faster access
    void configLocal(std::string path,std::string outPath = "");

    //Search the local data files and find the full gene stored in the local files
    std::string fullGene(std::string ncbiID, std::string dbName = "all");

    //Search the local files and if a more complete from is found then update the GeneData object list
    std::vector<GeneData> completeGenes(std::vector<GeneData> geneSet);

    //Fragment the given genes' gene sequence by index
    std::vector<GeneData> multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex);

    //Fragment the given genes' gene sequence by strings
    std::vector<GeneData> multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString);

    //Fragment the given gene's gene sequence by index
    GeneData indexFragment(GeneData inpData, int strtIndex, int endIndex);

    //Fragment the given gene's gene sequence by strings
    GeneData stringFragment(GeneData inpData, std::string strtString, std::string endString);

    //Train with an MSA image + gene sequences to identify particular letters
    void trainImageData(std::vector<std::string> genes, std::string path);

    //Delete all the trained data
    void removeImgTrainData();


    //******************************************************************************************************************
    //Complete Conversions
    //******************************************************************************************************************

    // Conversion of the Input Data to objects of standard format
    std::vector<std::string> objConvert(std::string inPath, std::string origin, std::string outType, int thrdCnt=10);

    //Main file conversion function
    void fileConvert(std::string inPath, std::string outPath, std::string origin, std::string outType, int thrdCnt=10);
}

#endif