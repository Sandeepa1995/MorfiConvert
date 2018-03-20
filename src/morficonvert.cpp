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

//Input Conversions
// Conversion of the Input Data to GeneData objects
std::vector<GeneData> morfiConvert::getInData(std::string path, std::string origin) {
    if (origin == "NCBI-BLAST") {
        return File2Obj::getInNCBIBLASTData(path);
    } else if (origin == "BLAST") {
        return File2Obj::getInCustomBLASTData(path);
    } else if (origin == "MSA") {
        return File2Obj::getInCustomMSAData(path);
    } else if (origin == "FASTA") {
        return File2Obj::getInFASTAData(path);
    } else if (origin == "FIND") {
        return morfiConvert::getInData(path, SupMorfi::identify(path));
    } else {
        //Throw exception if an unknown origin type
        throw UnknownInpException();
    }
}

//Output Conversions
// Conversion of GeneData objects to the requested format
std::vector<std::string> morfiConvert::giveOutData(std::vector<GeneData> res, std::string outType) {
    return Obj2Obj::giveOutData(res, outType);
}

void morfiConvert::writeOutData(std::vector<GeneData> res, std::string path, std::string outType){
    Obj2File::writeOutData(res,path,outType);
}

//Intermediaryy Functions
//Identify type
std::string morfiConvert::identify(std::string path){
    return SupMorfi::identify(path);
}

//Configure files for faster access
void morfiConvert::configLocal(std::string path){
    GeneComplete::configureFile(path);
}

std::string morfiConvert::completeGene(std::string ncbiID, std::string path){
    return GeneComplete::complete(ncbiID,path);
}

std::vector<GeneData> morfiConvert::completeGenes(std::vector<GeneData> geneSet){
    for (GeneData gene:geneSet){
        std::string apparentGene = GeneComplete::complete(SupMorfi::getIdentifier(gene.getDetails()),"all");
        if(apparentGene!="" && apparentGene.length()>gene.getGene().length()){
            gene.setGene(apparentGene);
        }
    }
    return geneSet;
}

std::vector<GeneData> morfiConvert::multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex){
    return GeneFragment::multiIndexFragment(inpData,strtIndex,endIndex);
}

std::vector<GeneData> morfiConvert::multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString){
    return GeneFragment::multiStringFragment(inpData,strtString,endString);
}

GeneData morfiConvert::indexFragment(GeneData inpData, int strtIndex, int endIndex){
    return GeneFragment::indexFragment(inpData, strtIndex, endIndex);
}

GeneData morfiConvert::stringFragment(GeneData inpData, std::string strtString, std::string endString){
    return GeneFragment::stringFragment(inpData, strtString, endString);
}


//Completed Functions
// Forced / Identified conversion of the Input Data to objects of standard format
std::vector<std::string> morfiConvert::objConvert(std::string inPath, std::string origin, std::string outType) {
    return Obj2Obj::giveOutData(morfiConvert::getInData(inPath, origin), outType);
}

//Main forced file conversion function
void morfiConvert::fileConvert(std::string inPath, std::string outPath, std::string origin, std::string outType) {
    Obj2File::writeOutData(morfiConvert::getInData(inPath, origin), outPath, outType);
}

////Main forced file conversion function
//void morfiConvert::MSAConvert(std::string inPath) {
//    File2Obj::getInMSADara(inPath);
//}

