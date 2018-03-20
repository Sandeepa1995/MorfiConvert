//
// Created by Damitha on 3/10/2018.
//

#include "obj2Obj.h"

#include <iostream>
#include <vector>

#include "GeneData.h"
#include "supplimentMorfi.h"

// Exception (super) class for MorfiConvert
class MorfiObj2ObjException : public std::exception {
public:
    virtual const char *what() const throw() {
        return "Error occured when converting from GeneData objects to the output set of objects";
    }
};

//Incorrect input type exception
class UnknownOutException : public MorfiObj2ObjException {
public:
    virtual const char *what() const throw() {
        return "Unsupported output data format type";
    }
};

//Return the data as a set of Vectors
std::vector<std::string> Obj2Obj::giveOutData(std::vector<GeneData> res, std::string outType = "FASTA") {
    //Initialize a output vector set
    std::vector<std::string> outVec;

#if defined(_OPENMP)
    //Check if output data type is FASTA
    if (outType == "FASTA") {
        //iterate through the GeneData set and write the data to the output vector set
#pragma  omp parallel for num_threads(10)
//        for (auto &value:res) {
        for (int x=0;x<res.size();x++) {
            GeneData value = res.at(x);
            //Temporary string to store the result that is being built
            std::string sTemp = "";
            sTemp += (">" + value.getDetails());
            //FASTA should be split into lines of a maximum of 80 characters
            int i = 0;
            for (char &c : value.getGene()) {
                if (i % 80 == 0) {
                    sTemp += "\n";
                }
                i += 1;
                sTemp += c;
            }
            outVec.push_back(sTemp);
        }
        return outVec;
        //Check if the output data type is the list of NCBI identifiers
        //!!! CONFIRM IF TO RETURN ONLY 1 OR ALL !!!
    } else if (outType == "NCBI-IDs") {
        //Iterate through the list of gene sequences
#pragma  omp parallel for num_threads(10)
//        for (auto &value:res) {
        for (int x=0;x<res.size();x++) {
            GeneData value = res.at(x);
            //Temporary variable to hold the gene description
            std::string gene = value.getDetails();
            outVec.push_back(SupMorfi::getIdentifier(gene));
//            //Iterate through the gene sequence's details to search for a '.' which is an essential component of the NCBI ID
//            for (int i = 0; i < gene.length(); ++i) {
//                if (gene[i] == '.') {
//                    //Variables to hold the starting and ending positions
//                    int pre = i - 1;
//                    int post = i + 1;
//                    //Iterate forward till a white space (Number / Version portion of the identifier)
//                    while (gene[post] != ' ' && post <= gene.length()) {
//                        post += 1;
//                    }
//                    //Iterate backward till white space (Actual ID)
//                    while (gene[pre] != ' ' && pre >= 0) {
//                        pre -= 1;
//                    }
//                    //Check if the selected region is matching with the NCBI ID format and then add to the output vector set
//                    if (isdigit(gene[post - 1]) && isalnum(gene[pre + 1])) {
//                        outVec.push_back(gene.substr(pre + 1, post - pre));
//                    }
//                }
//            }
        }

        return outVec;
    } else {
        //Throw error for unknown type
        throw UnknownOutException();
    }
#else
    if (outType == "FASTA") {
        //iterate through the GeneData set and write the data to the output vector set
        for (auto &value:res) {
            //Temporary string to store the result that is being built
            std::string sTemp = "";
            sTemp += (">" + value.getDetails());
            //FASTA should be split into lines of a maximum of 80 characters
            int i = 0;
            for (char &c : value.getGene()) {
                if (i % 80 == 0) {
                    sTemp += "\n";
                }
                i += 1;
                sTemp += c;
            }
            outVec.push_back(sTemp);
        }
        return outVec;
        //Check if the output data type is the list of NCBI identifiers
        //!!! CONFIRM IF TO RETURN ONLY 1 OR ALL !!!
    } else if (outType == "NCBI-IDs") {
        //Iterate through the list of gene sequences
        for (auto &value:res) {
            //Temporary variable to hold the gene description
            std::string gene = value.getDetails();
            //Iterate through the gene sequence's details to search for a '.' which is an essential component of the NCBI ID
            for (int i = 0; i < gene.length(); ++i) {
                if (gene[i] == '.') {
                    //Variables to hold the starting and ending positions
                    int pre = i - 1;
                    int post = i + 1;
                    //Iterate forward till a white space (Number / Version portion of the identifier)
                    while (gene[post] != ' ' && post <= gene.length()) {
                        post += 1;
                    }
                    //Iterate backward till white space (Actual ID)
                    while (gene[pre] != ' ' && pre >= 0) {
                        pre -= 1;
                    }
                    //Check if the selected region is matching with the NCBI ID format and then add to the output vector set
                    if (isdigit(gene[post - 1]) && isalnum(gene[pre + 1])) {
                        outVec.push_back(gene.substr(pre + 1, post - pre));
                    }
                }
            }
        }

        return outVec;
    } else {
        //Throw error for unknown type
        throw UnknownOutException();
    }
#endif
}