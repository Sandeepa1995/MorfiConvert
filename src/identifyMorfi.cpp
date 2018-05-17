//
// Created by Damitha on 5/12/2018.
//

#include "identifyMorfi.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include "include/json/json/json.h"
#include "GeneData.h"
#include "supplimentMorfi.h"

// Exception (super) class for MorfiConvert File2Obj
class MorfiSupException : public std::exception {
public:
    virtual const char *what() const throw() {
        return "Error occurred in supplementary function.";
    }
};

//Input file opening custom exception
class FileInpException : public MorfiSupException {
public:
    virtual const char *what() const throw() {
        return "Error in opening input file. Please check if the given path is correct.";
    }
};

//Incorrect input file type exception
class UnknownInpFileException : public MorfiSupException {
public:
    virtual const char *what() const throw() {
        return "Unsupported input file format";
    }
};

//Check if it is of type NCBI BLAST ( if it is a .txt file)
bool isNCBIBLAST_TXT(std::string path) {

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        bool isNCBIBLAST = false;
        //String to hold every file-line read
        std::string s;

        //Iterate through the file
        while (getline(infile, s)) {
            //Once "BLAST" keyword is found iterate from that point on to find matches
            if (s.length() > 0 && s.substr(0, 5) == "BLAST") {
                while (getline(infile, s)) {
                    //Locate FASTA description
                    if (s.length() > 0 && s[0] == '>') {
                        while (getline(infile, s)) {
                            //Locate matched gene sequence line
                            if (s.length() > 0 && s.substr(0, 5) == "Sbjct") {
                                //String stream to iterate through a line and find the necessary sections
                                std::istringstream buf(s);
                                std::istream_iterator<std::string> beg(buf), end;

                                std::vector<std::string> tokens(beg, end);
                                //Check if the 3rd element in the line matches a gene sequence pattern
                                if (tokens[2].length() > 40) {
                                    //Identify as NCBI BLAST
                                    isNCBIBLAST = true;
                                    break;
                                }
                            }
                        }
                    }

                    //Break if identified as a NCBI BLAST file.
                    if (isNCBIBLAST) {
                        break;
                    }
                }
            }
            //Break if identified as a NCBI BLAST file.
            if (isNCBIBLAST) {
                break;
            }
        }

        //Return the result
        return isNCBIBLAST;

    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//Check if it is of type NCBI BLAST ( if it is a .json file)
bool isNCBIBLAST_JSON(std::string path) {
    //Variable to hold all json values
    Json::Value root;

    //Get file and parse the values in it as json and assign to the holder variable
    std::ifstream infile(path);
    infile >> root;

    //Check if the input format is NCBI BLAST
    bool isNCBIBLAST = false;

    //Check each subsection of the json file to see if it matches the format of the NCBI BLAST result in json
    if (root.isMember("BlastOutput2")) {
        if (root["BlastOutput2"].isMember("report")) {
            if (root["BlastOutput2"]["report"].isMember("results")) {
                if (root["BlastOutput2"]["report"]["results"].isMember("search")) {
                    if (root["BlastOutput2"]["report"]["results"]["search"].isMember("hits")) {
                        if (root["BlastOutput2"]["report"]["results"]["search"]["hits"][0].isMember("hsps")) {
                            if (root["BlastOutput2"]["report"]["results"]["search"]["hits"][0]["hsps"][0].isMember(
                                    "hseq")) {
                                std::string testGene = root["BlastOutput2"]["report"]["results"]["search"]["hits"][0]["hsps"][0]["hseq"].asString();
                                isNCBIBLAST = true;
                                //Check if there is any space in the FASTA gene sequence
                                for (char &c : testGene) {
                                    if (c == ' ') {
                                        isNCBIBLAST = false;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (!isNCBIBLAST) {
            if (root["BlastOutput2"][0].isMember("report")) {
                if (root["BlastOutput2"][0]["report"].isMember("results")) {
                    if (root["BlastOutput2"][0]["report"]["results"].isMember("search")) {
                        if (root["BlastOutput2"][0]["report"]["results"]["search"].isMember("hits")) {
                            if (root["BlastOutput2"][0]["report"]["results"]["search"]["hits"][0].isMember("hsps")) {
                                if (root["BlastOutput2"][0]["report"]["results"]["search"]["hits"][0]["hsps"][0].isMember(
                                        "hseq")) {
                                    std::string testGene = root["BlastOutput2"][0]["report"]["results"]["search"]["hits"][0]["hsps"][0]["hseq"].asString();
                                    isNCBIBLAST = true;
                                    //Check if there is any space in the FASTA gene sequence
                                    for (char &c : testGene) {
                                        if (c == ' ') {
                                            isNCBIBLAST = false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return isNCBIBLAST;
}

//Check if it is of type NCBI BLAST ( if it is a .xml file)
bool isNCBIBLAST_XML(std::string path) {
    std::ifstream infile;

    //Check if the input format is NCBI BLAST
    bool isNCBIBLAST = false;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        //String to hold every file-line read
        std::string s;

        //Read entire file, ( FASTA section by FASTA section )
        while (getline(infile, s)) {
            if (s.substr(s.find_first_not_of(" "), s.find_last_of(" ")) == "<BlastOutput2") {
                while (getline(infile, s)) {
                    if (s.substr(s.find_first_not_of(" "), s.find_last_of(" ")) == "<Hit>") {
                        while (getline(infile, s)) {
                            if (s.substr(s.find_first_not_of(" "), 4) == "<id>") {
                                while (getline(infile, s)) {
                                    if (s.substr(s.find_first_not_of(" "), 7) == "<title>") {
                                        while (getline(infile, s)) {
                                            if (s.substr(s.find_first_not_of(" "), 6) == "<hseq>") {
                                                std::string testGene = s.substr(s.find_first_of(">") + 1,
                                                                                s.find_last_of("<") -
                                                                                s.find_first_of(">") - 1);
                                                if (testGene.length() > 40) {
                                                    isNCBIBLAST = true;
                                                    break;
                                                }
                                                //Check if there is any space in the FASTA gene sequence
                                                for (char &c : testGene) {
                                                    if (c == ' ') {
                                                        isNCBIBLAST = false;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if (isNCBIBLAST) {
                                        break;
                                    }
                                }
                            }
                            if (isNCBIBLAST) {
                                break;
                            }
                        }
                    }
                    if (isNCBIBLAST) {
                        break;
                    }
                }
            }
            if (isNCBIBLAST) {
                break;
            }
        }
        infile.clear();
        infile.close();
    } else {
        throw FileInpException();
    }

    //If there was no error in opening the file
    if (infile.good()) {
        //String to hold every file-line read
        std::string s;

        //Read entire file, ( FASTA section by FASTA section )
        while (getline(infile, s)) {
            if (s.substr(s.find_first_not_of(" "), s.find_last_of(" ")) == "<BlastOutput2") {
                while (getline(infile, s)) {
                    if (s.substr(s.find_first_not_of(" "), s.find_last_of(" ")) == "<Hit>") {
                        while (getline(infile, s)) {
                            //2 values are checked due to one being in the protien BLAST while the other is in the
                            // nucleic acid BLAST
                            if (s.substr(s.find_first_not_of(" "), 4) == "<id>" ||
                                s.substr(s.find_first_not_of(" "), 8) == "<Hit_id>") {
                                while (getline(infile, s)) {
                                    //2 values are checked due to one being in the protien BLAST while the other is in
                                    // the nucleic acid BLAST
                                    if (s.substr(s.find_first_not_of(" "), 6) == "<hseq>" ||
                                        s.substr(s.find_first_not_of(" "), 10) == "<Hsp_hseq>") {
                                        std::string testGene = s.substr(s.find_first_of(">") + 1,
                                                                        s.find_last_of("<") -
                                                                        s.find_first_of(">") - 1);
                                        if (testGene.length() > 40) {
                                            isNCBIBLAST = true;
                                            break;
                                        }
                                        //Check if there is any space in the FASTA gene sequence
                                        for (char &c : testGene) {
                                            if (c == ' ') {
                                                isNCBIBLAST = false;
                                            }
                                        }
                                    }
                                    if (isNCBIBLAST) {
                                        break;
                                    }
                                }
                            }
                            if (isNCBIBLAST) {
                                break;
                            }
                        }
                    }
                    if (isNCBIBLAST) {
                        break;
                    }
                }
            }
            if (isNCBIBLAST) {
                break;
            }
        }
        infile.clear();
        infile.close();
    } else {
        throw FileInpException();
    }

    return isNCBIBLAST;
}

//Check if it is of type GCG
bool isGCG(std::string path) {

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        bool isGCG = false;
        //String to hold every file-line read
        std::string s;

        //Iterate through the file
        while (getline(infile, s)) {
            if (s.length() > 0 && s.substr(0, 13) == "!!NA_SEQUENCE" ||
                s.length() > 0 && s.substr(0, 13) == "!!AA_SEQUENCE") {
                //Identify as GCG
                isGCG = true;
                break;
            }
        }
        //Return the result
        return isGCG;
    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//Check if it is of type FASTA
bool isFASTA(std::string path) {

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        bool isFASTA = false;
        //String to hold every file-line read
        std::string s;

        //Iterate through the file
        while (getline(infile, s)) {
            if (s.length() == 0 || isFASTA) {
                break;
            }
            if (s.length() > 0 && s.substr(0, 1) == ">") {
                while (getline(infile, s)) {
                    if (s.length() == 0) {
                        break;
                    }
                    if (s.length() > 0 && s.substr(0, 1) == ">") {
                        //Identify as FASTA
                        isFASTA = true;
                        break;
                    }
                }
            }
        }
        //Return the result
        return isFASTA;
    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//Check if it is of type EMBL
bool isEMBL(std::string path) {

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        bool isEMBL = false;
        //String to hold every file-line read
        std::string s;

        //Iterate through the file
        while (getline(infile, s)) {
            if (s.length() > 0 && s.substr(0, 2) == "ID") {
                while (getline(infile, s)) {
                    if (s.length() > 0 && s.substr(0, 2) == "XX") {
                        while (getline(infile, s)) {
                            if (s.length() > 0 && s.substr(0, 2) == "DE") {
                                while (getline(infile, s)) {
                                    if (s.length() > 0 && s.substr(0, 2) == "XX") {
                                        while (getline(infile, s)) {
                                            if (s.length() > 0 && s.substr(0, 2) == "SQ") {
                                                //Identify as EMBL
                                                isEMBL = true;
                                                break;
                                            }
                                        }
                                    }
                                    if(isEMBL){
                                        break;
                                    }
                                }
                            }
                            if(isEMBL){
                                break;
                            }
                        }
                    }
                    if(isEMBL){
                        break;
                    }
                }
            }
            if(isEMBL){
                break;
            }
        }

        //Return the result
        return isEMBL;
    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//Check if it is of type GenBank
bool isGenBank(std::string path) {

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        bool isGenBank = false;
        //String to hold every file-line read
        std::string s;

        //Iterate through the file
        while (getline(infile, s)) {
            if (s.length() > 0 && s.substr(0, 5) == "LOCUS") {
                while (getline(infile, s)) {
                    if (s.length() > 0 && s.substr(0, 10) == "DEFINITION") {
                        while (getline(infile, s)) {
                            if (s.length() > 0 && s.substr(0, 6) == "ORIGIN") {
                                while (getline(infile, s)) {
                                    if (s.length() > 0 && s.substr(0, 2) == "//") {
                                        //Identify as GenBank
                                        isGenBank = true;
                                        break;
                                    }
                                }
                            }
                            if(isGenBank){
                                break;
                            }
                        }
                    }
                    if(isGenBank){
                        break;
                    }
                }
            }
            if(isGenBank){
                break;
            }
        }

        //Return the result
        return isGenBank;
    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//Check if it is of type PIR / NBRF
bool isPIRorNBRF(std::string path) {

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        bool isPIR = false;
        //String to hold every file-line read
        std::string s;

        //Iterate through the file
        while (getline(infile, s)) {
            if (s.length() > 0 && s.substr(0, 4) == ">P1;" ||
                s.length() > 0 && s.substr(0, 4) == ">D1;") {
                //Identify as PIR
                isPIR = true;
                break;
            }
        }
        //Return the result
        return isPIR;
    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//Check if it is of type PHYLIP
bool isPHYLIP(std::string path) {

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        bool isPHYLIP = false;
        //String to hold every file-line read
        std::string s;

        getline(infile, s);
        if(isdigit(s[0])){
            getline(infile, s);
            if(s[0]!=' '){
                getline(infile, s);
                if(s[0] == ' '){
                    //Identify as PHYLIP
                    isPHYLIP = true;
                }
            }
        }
        //Return the result
        return isPHYLIP;
    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//Check if it is of type CLUSTAL
bool isCLUSTAL(std::string path) {

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        bool isCLUSTAL = false;
        //String to hold every file-line read
        std::string s;

        getline(infile, s);
        if(s.length() > 0 && s.substr(0, 7) == "CLUSTAL"){
            //Identify as CLUSTAL
            isCLUSTAL = true;
        }
        //Return the result
        return isCLUSTAL;
    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//Check if it is of type FASTA Report
bool isFASTARep(std::string path) {

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        bool isFASTARep = false;
        //String to hold every file-line read
        std::string s;

        //Check first line for "FASTA"
        getline(infile, s);
        if(s.length() > 0 && s.substr(0, 5) == "FASTA"){
            //Identify as CLUSTAL
            isFASTARep = true;
        }

        //Check second line for "FASTA" as well
        getline(infile, s);
        if(s.length() > 0 && s.substr(0, 5) == "FASTA"){
            //Identify as CLUSTAL
            isFASTARep = true;
        }
        //Return the result
        return isFASTARep;
    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//Identify the data type of a given file
std::string IdentifyMorfi::identify(std::string path) {
    //Check the file content 1 by 1 to see if it matches with standard formats

    //GCG
    if (isGCG(path)) {
        return "GCG";
    }

    //FASTA
    if (isFASTA(path)) {
        return "FASTA";
    }

    //EMBL
    if (isEMBL(path)) {
        return "EMBL";
    }

    //GenBank
    if (isGenBank(path)) {
        return "GenBank";
    }

    //PIR/NBRF
    if (isPIRorNBRF(path)) {
        return "PIR";
    }

    //PHYLIP
    if (isPHYLIP(path)) {
        return "PHYLIP";
    }

    //CLUSTAL
    if (isCLUSTAL(path)) {
        return "Clustral";
    }

    //FASTA-Report
    if (isFASTARep(path)) {
        return "FASTA-Report";
    }

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

    //NCBI-BLAST
    if (ext == "txt") {
        if (isNCBIBLAST_TXT(path)) {
            return "NCBI-BLAST";
        }

    } else if (ext == "json") {
        if (isNCBIBLAST_JSON(path)) {
            return "NCBI-BLAST";
        }

    } else if (ext == "xml") {
        if (isNCBIBLAST_XML(path)) {
            return "NCBI-BLAST";
        }
    }
    return "Unknown";
}