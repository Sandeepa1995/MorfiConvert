//
// Created by Damitha on 3/10/2018.
//

#include "supplimentMorfi.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include "include/json/json/json.h"

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

std::string SupMorfi::getExtension(std::string path) {
    std::string lastDash = path.substr(path.find_last_of("\\") + 1);
    return lastDash.substr(lastDash.find_last_of(".") + 1);
}

bool SupMorfi::isAlpha(std::string s) {
    for (char c:s) {
        int cha = (int) c;
        if (!((cha <= (int) 'z' && cha >= (int) 'a') || (cha <= (int) 'Z' && cha >= (int) 'A'))) {
            return false;
        }
    }
    return true;
}

std::string SupMorfi::getIdentifier(std::string s) {
    for (int i = 0; i < s.length(); ++i) {
        if (s[i] == '.') {
            //Variables to hold the starting and ending positions
            int pre = i - 1;
            int post = i + 1;
            //Iterate forward till a white space (Number / Version portion of the identifier)
            while (s[post] != ' ' && post <= s.length()) {
                post += 1;
            }
            //Iterate backward till white space (Actual ID)
            while (s[pre] != ' ' && pre >= 0) {
                pre -= 1;
            }
            //Check if the selected region is matching with the NCBI ID format and then add to the output vector set
            if (isdigit(s[post - 1]) && isalnum(s[pre + 1])) {
                return (s.substr(pre + 1, post - pre));
            }
        }
    }
}

std::string SupMorfi::identify(std::string path) {
    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

    //txt identification
    if (ext == "txt") {
        //Check if it is of type NCBI BLAST
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

                        if (isNCBIBLAST) {
                            break;
                        }
                    }
                }
                if (isNCBIBLAST) {
                    break;
                }
            }

            if (isNCBIBLAST) {
                return "NCBI-BLAST";
            }

            return "Unknwon";
        } else {
            //Throw error in file open
            throw FileInpException();
        }
    } else if (ext == "json") {
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
        }

        if (isNCBIBLAST) {
            return "NCBI-BLAST";
        } else {
            return "Unknown";
        }

    } else if (ext == "xml") {
        std::ifstream infile;

        //Open file
        infile.open(path);

        //If there was no error in opening the file
        if (infile.good()) {
            //String to hold every file-line read
            std::string s;

            //Check if the input format is NCBI BLAST
            bool isNCBIBLAST = false;

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

            if (isNCBIBLAST) {
                return "NCBI-BLAST";
            } else {
                return "Unknown";
            }
        } else {
            throw FileInpException();
        }
    } else {
        throw UnknownInpFileException();
    }
}