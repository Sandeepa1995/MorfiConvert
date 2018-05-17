#include "morficonvert.h"

#include <iostream>
#include <vector>

#include "GeneData.h"
#include "file2Obj.h"
#include "obj2Obj.h"
#include "obj2File.h"
#include "supplimentMorfi.h"
#include "geneComplete.h"
#include "geneFragment.h"
#include "image2Text.h"
#include "identifyMorfi.h"

//Exceptions
// Exception (super) class for MorfiConvert
class MorfiConvertException : public std::exception {
public:
    virtual const char *what() const throw() {
        return "Unknown error occurred.";
    }
};

//Incorrect input type exception
class UnknownInpException : public MorfiConvertException {
public:
    virtual const char *what() const throw() {
        return "Unsupported input data format type";
    }
};

//**********************************************************************************************************************
//Input Conversions
//**********************************************************************************************************************

// Conversion of the Input Data to GeneData objects
std::vector<GeneData> morfiConvert::getInData(std::string path, std::string origin, int thrdCnt) {
    if (origin == "NCBI-BLAST") {
        return File2Obj::getInNCBIBLASTData(path, thrdCnt);
    } else if (origin == "BLAST") {
        return File2Obj::getInCustomBLASTData(path, thrdCnt);
    } else if (origin == "MSA") {
        return File2Obj::getInCustomMSAData(path);
    } else if (origin == "FASTA") {
        return File2Obj::getInFASTAData(path, thrdCnt);
    } else if (origin == "GCG") {
        return File2Obj::getInGCGData(path, thrdCnt);
    } else if (origin == "EMBL") {
        return File2Obj::getInEMBLData(path, thrdCnt);
    } else if (origin == "GenBank") {
        return File2Obj::getInGenebankData(path, thrdCnt);
    } else if (origin == "PIR" || origin == "NBRF") {
        return File2Obj::getInPIRData(path, thrdCnt);
    } else if (origin == "PHYLIP") {
        return File2Obj::getInPHYLIPData(path, thrdCnt);
    } else if (origin == "Clustral") {
        return File2Obj::getInClustralData(path);
    } else if (origin == "MSF") {
        return File2Obj::getInMSFData(path);
    } else if (origin == "FASTA-Report") {
        return File2Obj::getInFASTARepData(path, thrdCnt);
    } else if (origin == "FIND") {
        //If the input type is given to be "FIND" then identify the given type via the type identification function
        return morfiConvert::getInData(path, IdentifyMorfi::identify(path), thrdCnt);
    } else {
        //Throw exception if an unknown origin type
        throw UnknownInpException();
    }
}

// Conversion of the Input MSA Data to GeneData objects (Image files)
std::vector<GeneData> morfiConvert::getInImageData(std::vector<std::string> genes, std::string path){
    return Image2Text::getInImgMSAData(genes,path);
}

//**********************************************************************************************************************
//Output Conversions
//**********************************************************************************************************************

// Conversion of GeneData objects to the requested string format
std::vector<std::string> morfiConvert::giveOutData(std::vector<GeneData> res, std::string outType, int thrdCnt) {
    return Obj2Obj::giveOutData(res, outType, thrdCnt);
}

//Write the given GeneData objects to a specific file format
void morfiConvert::writeOutData(std::vector<GeneData> res, std::string path, std::string outType, int thrdCnt) {
    //PHYLIP files would contain the additional "<Number of genes> <Longest gene sequence>" line
    if (outType == "PHYLIP") {
        //Convert the internal data
        std::vector<std::string> inpStrVec = morfiConvert::giveOutData(res, outType, thrdCnt);

        //Compute and add the first line
        inpStrVec.insert(inpStrVec.begin(),
                         (std::to_string(res.size()) + " " + std::to_string(SupMorfi::getLongestGeneLength(res))));
        //Write the data
        Obj2File::writeOutData(inpStrVec, path, outType);

    } else {
        //Convert to the needed type then write the content
        Obj2File::writeOutData(morfiConvert::giveOutData(res, outType, thrdCnt), path, outType);
    }

}
//**********************************************************************************************************************
//Intermediary / supplimentary Functions
//**********************************************************************************************************************
//Identify type
std::string morfiConvert::identify(std::string path) {
    return IdentifyMorfi::identify(path);
}

//Configure files for faster access
void morfiConvert::configLocal(std::string path, std::string outPath) {
    GeneComplete::configureFile(path,outPath);
}

//Search the local data files and find the full gene stored in the local files
std::string morfiConvert::fullGene(std::string ncbiID, std::string dbName) {
    return GeneComplete::complete(ncbiID, dbName);
}

//Search the local files and if a more complete from is found then update the GeneData object list
std::vector<GeneData> morfiConvert::completeGenes(std::vector<GeneData> geneSet) {
    for (GeneData gene:geneSet) {
        std::string apparentGene = GeneComplete::complete(SupMorfi::getIdentifier(gene.getDetails()), "all");
        //Only update the gene sequence if the gene found in the local storage is longer than the initial gene
        if (apparentGene != "" && apparentGene.length() > gene.getGene().length()) {
            gene.setGene(apparentGene);
        }
    }
    return geneSet;
}

//Fragment the given genes' gene sequence by index
std::vector<GeneData> morfiConvert::multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex) {
    return GeneFragment::multiIndexFragment(inpData, strtIndex, endIndex);
}

//Fragment the given genes' gene sequence by strings
std::vector<GeneData>
morfiConvert::multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString) {
    return GeneFragment::multiStringFragment(inpData, strtString, endString);
}

//Fragment the given gene's gene sequence by index
GeneData morfiConvert::indexFragment(GeneData inpData, int strtIndex, int endIndex) {
    return GeneFragment::indexFragment(inpData, strtIndex, endIndex);
}

//Fragment the given gene's gene sequence by strings
GeneData morfiConvert::stringFragment(GeneData inpData, std::string strtString, std::string endString) {
    return GeneFragment::stringFragment(inpData, strtString, endString);
}

//Train with an MSA image + gene sequences to identify particular letters
void morfiConvert::trainImageData(std::vector<std::string> genes, std::string filePath){
    Image2Text::trainImgMSAData(genes,filePath);
}

//Delete all the trained data
void morfiConvert::removeImgTrainData(){
    Image2Text::removeTrainData();
}

//**********************************************************************************************************************
//Complete Conversions
//**********************************************************************************************************************

// Forced / Identified conversion of the Input Data to objects of standard format
std::vector<std::string>
morfiConvert::objConvert(std::string inPath, std::string origin, std::string outType, int thrdCnt) {
    return Obj2Obj::giveOutData(morfiConvert::getInData(inPath, origin, thrdCnt), outType, thrdCnt);
}

//Main forced file conversion function
void morfiConvert::fileConvert(std::string inPath, std::string outPath, std::string origin, std::string outType,
                               int thrdCnt) {
    if (outType == "PHYLIP") {
        std::vector<GeneData> res = morfiConvert::getInData(inPath, origin, thrdCnt);
        std::vector<std::string> inpStrVec = morfiConvert::giveOutData(res, outType, thrdCnt);
        inpStrVec.insert(inpStrVec.begin(),
                         (std::to_string(res.size()) + " " + std::to_string(SupMorfi::getLongestGeneLength(res))));
        Obj2File::writeOutData(inpStrVec, outPath, outType);
    } else {
        Obj2File::writeOutData(
                morfiConvert::giveOutData(morfiConvert::getInData(inPath, origin, thrdCnt), outType, thrdCnt), outPath,
                outType);
    }
}

