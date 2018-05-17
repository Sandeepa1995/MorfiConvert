#include "file2Obj.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <memory>
#include <omp.h>

#include "include/json/json/json.h"
#include "GeneData.h"
#include "supplimentMorfi.h"

//Custom exceptions

// Exception (super) class for MorfiConvert File2Obj
class MorfiFile2ObjException : public std::exception {
public:
    virtual const char *what() const throw() {
        return "Error when converting the file to a set of GeneData objects.";
    }
};

//Input file opening custom exception
class FileInpException : public MorfiFile2ObjException {
public:
    virtual const char *what() const throw() {
        return "Error in opening input file. Please check if the given path is correct.";
    }
};

//Unknown input file type exception
class UnknownInpFileException : public MorfiFile2ObjException {
public:
    virtual const char *what() const throw() {
        return "Unsupported input file format";
    }
};

//Incorrect input file type exception
class IncorrectInpDataException : public MorfiFile2ObjException {
public:
    virtual const char *what() const throw() {
        return "Incorrect input file data format";
    }
};

//Function to read the file in the given path and return it as a vector of strings
std::vector<std::string> getInpFileAsStrings(std::string path) {
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //String to hold every file-line read
    std::vector<std::string> sVec;

    //If there was no error in opening the file
    if (infile.good()) {

        //Temporary storage string
        std::string sTemp;

        // Add all the content in the file to sVec (To incorporate parallelism when processing the data)
        while (getline(infile, sTemp)) {
            sVec.push_back(sTemp);
        }

        //Close and clear input file stream once done
        infile.close();
        infile.clear();
    } else {
        //Throw error in file open
        throw FileInpException();
    }

    //Return the result
    return sVec;
}

//Most of the following conversions will have the path of the data file as well as the threads that the conversion
// process is to run on. If 1 then will be run in sequence.

//The NCBI BLAST data file conversion function. Supports .txt, .json and .xml formats.
std::vector<GeneData> File2Obj::getInNCBIBLASTData(std::string path, int thrdCnt) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

    //txt extraction
    if (ext == "txt") {

        //Get the input file as a vector string to support parallel processing.
        std::vector<std::string> sVec = getInpFileAsStrings(path);

        //OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
        for (int i = 0; i < sVec.size(); i++) {
            //Check FASTA description start
            std::string s1 = sVec.at(i);
            if (s1.length() > 0 && s1[0] == '>') {
                int pointerFASTA = 1;
                //New GeneData Object
                GeneData newInData("", "");
                //Add the current line to the description
                newInData.setDetails(s1.substr(1));

                int stopJ;
                //Iterate and add to description till 'Length' field
                for (int j = i + 1; j < sVec.size(); j++) {
                    std::string s2 = sVec.at(j);
                    if (s2.length() > 0 && s2[0] == 'L' && s2[1] == 'e' && s2[2] == 'n') {
                        stopJ = j;
                        break;
                    } else {
                        newInData.setDetails(newInData.getDetails() + s2);
                    }
                }

                //Iterate through the FASTA sequences
                for (int k = stopJ + 1; k < sVec.size(); k++) {
                    std::string s3 = sVec.at(k);
                    //Check if subject FASTA
                    if (s3.length() > 0 && s3[0] == 'S' && s3[1] == 'b' && s3[2] == 'j' && s3[3] == 'c' &&
                        s3[4] == 't') {

                        //String stream to iterate through a line and find the necessary sections
                        std::istringstream buf(s3);
                        std::istream_iterator<std::string> beg(buf), end;

                        std::vector<std::string> tokens(beg, end);
                        std::stringstream toInt(tokens[1]);
                        int currentPointer = 1;
                        toInt >> currentPointer;
                        //Add to the FASTA sequence
                        //Check pointer in FASTA, in case a description is missing.
                        if (currentPointer >= pointerFASTA) {
                            newInData.setGene(newInData.getGene() + tokens[2]);
                        } else {
                            if (newInData.getGene().length() > 0) {
#pragma omp critical
                                {
                                    inputData.push_back(newInData);
                                }
                            }
                            break;
                        }
                        pointerFASTA = currentPointer;

                    } else if (s3.length() > 0 && ((s3[0] == '>') ||
                                                   (s3[0] == ' ' && s3[1] == ' ' && s3[2] == 'D' &&
                                                    s3[3] == 'a' &&
                                                    s3[4] == 't'))) {
                        //Add the GeneData to the inputData set once a new FASTA description is discovered and recognized
                        if (newInData.getGene().length() > 0) {
                            //Make adding the GeneData object to the output vector thread safe
#pragma omp critical
                            {
                                inputData.push_back(newInData);
                            }
                        }
                        break;
                    }
                }
            }
        }
        //Return the transformed data
        if (inputData.size() > 0) {
            return inputData;
        } else {
            throw IncorrectInpDataException();
        }
    } else if (ext == "json") {
        //Variable to hold all json values
        Json::Value root;

        //Get file and parse the values in it as json and assign to the holder variable
        std::ifstream infile(path);
        infile >> root;

        Json::Value res = "";
        //Array that the BLAST results iterate over
        if (root.isMember("BlastOutput2")) {
            //Two separate situations considering protein and nucleic acid BLAST results.
            //Protein
            if (root["BlastOutput2"].isMember("report")) {
                if (root["BlastOutput2"]["report"].isMember("results")) {
                    if (root["BlastOutput2"]["report"]["results"].isMember("search")) {
                        if (root["BlastOutput2"]["report"]["results"]["search"].isMember("hits")) {
                            res = root["BlastOutput2"]["report"]["results"]["search"]["hits"];
                        }
                    }
                }
                //Nucleic acid
            } else if (root["BlastOutput2"][0].isMember("report")) {
                if (root["BlastOutput2"][0]["report"].isMember("results")) {
                    if (root["BlastOutput2"][0]["report"]["results"].isMember("search")) {
                        if (root["BlastOutput2"][0]["report"]["results"]["search"].isMember("hits")) {
                            res = root["BlastOutput2"][0]["report"]["results"]["search"]["hits"];
                        }
                    }
                }
            }
        }

        //Iterate over BLAST results
        for (Json::ValueConstIterator it = res.begin(); it != res.end(); ++it) {
            //New GeneData object
            GeneData newInData("", "");

            //Variable to hold the description of the gene sequences in the BLAST file
            std::string geneDisc = "";

            //Array that the description iterates over *Some can have multiple*
            const Json::Value &geneSet = *it;

            //Iterate over descriptions
            for (Json::ValueConstIterator itlv2 = geneSet["description"].begin();
                 itlv2 != geneSet["description"].end(); ++itlv2) {
                //Object to get the description of a gene sequence
                const Json::Value &gene = *itlv2;
                //Extract the NCBI ID of the gene
                std::string exrt = gene["id"].asString().substr(0, gene["id"].asString().substr(0,
                                                                                                gene["id"].asString().find_last_of(
                                                                                                        "|")).length());
                std::string tempDetail = (exrt.substr(exrt.find_last_of('|') + 1) + " ");
                geneDisc += tempDetail;
                //Add the title to the gene sequence description
                geneDisc += (gene["title"].asString() + " ");
            }
            //Remove tailing whitespace
            newInData.setDetails(geneDisc.substr(0, geneDisc.length() - 1));
            //Set the FASTA sequence
            newInData.setGene(geneSet["hsps"][0]["hseq"].asString());
            //Add the gene data to the vector set
            if (SupMorfi::getIdentifier(newInData.getDetails()).length() > 0 && newInData.getGene().length() > 0) {
                inputData.push_back(newInData);
            }
        }

        //Return the transformed data
        if (inputData.size() > 0) {
            return inputData;
        } else {
            throw IncorrectInpDataException();
        }
    } else if (ext == "xml") {
        //String to hold every file-line read
        std::vector<std::string> sVec = getInpFileAsStrings(path);

//OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
        for (int i = 0; i < sVec.size(); i++) {
            std::string s1 = sVec.at(i);
            //Start of a gene sequence used for BLAST
            if (s1.substr(s1.find_first_not_of(" "), s1.find_last_of(" ")) == "<Hit>") {
                //New GeneData object
                GeneData newInData("", "");
                //To hold the description of the gene sequence being processed
                std::string genDisc = "";

                //iterate through all the gene sequence descriptions and FASTA sequence
                for (int j = i + 1; j < sVec.size(); j++) {
                    std::string s2 = sVec.at(j);
                    //If the FASTA sequence is located then add to GeneData and break from loop
                    if ((s2.substr(s2.find_first_not_of(" "), 6) == "<hseq>")) {
                        //Get substring to avoid taking xml tags
                        newInData.setGene(
                                s2.substr(s2.find_first_not_of(" ") + 6, s2.substr(s2.find_first_not_of(" ")).
                                        length() - 13));
                        break;
                    } else if ((s2.substr(s2.find_first_not_of(" "), 10) == "<Hsp_hseq>")) {
                        //Get substring to avoid taking xml tags
                        newInData.setGene(
                                s2.substr(s2.find_first_not_of(" ") + 10, s2.substr(s2.find_first_not_of(" ")).
                                        length() - 21));
                        break;
                    }
                    if (s2.substr(s2.find_first_not_of(" "), 6) == "</Hsp>") {
                        break;
                    }


                    //A description of a gene sequence, iterate within to get the NCBI identifier and the additional
                    //description. AVOID LEADING WHITE SPACES / REMEMBER TO TRIM THEM!!
                    if (s2.substr(s2.find_first_not_of(" "), s2.find_last_of(" ")) == "<HitDescr>") {
                        //The extraction of data from protien BLAST results
                        for (int k = j + 1; k < sVec.size(); k++) {
                            std::string s3 = sVec.at(k);
                            //Check <id> tag to get the identifier
                            if (s3.substr(s3.find_first_not_of(" "), 4) == "<id>") {
                                //To avoid taking tags
                                std::string singleGeneDesc = s3.substr(s3.find_first_not_of(" ") + 4,
                                                                       (s3.substr(s3.find_first_not_of(" "),
                                                                                  s3.find_last_not_of(
                                                                                          " "))).length() -
                                                                       9);
                                //Add the description to the GeneData object
                                singleGeneDesc = singleGeneDesc.substr(0, singleGeneDesc.find_last_of("|"));
                                newInData.setDetails(newInData.getDetails() +
                                                     singleGeneDesc.substr(singleGeneDesc.find_last_of("|") + 1) +
                                                     " ");
                            }

                            //Check <title> tag for the additional identifier. (The scientific name etc.)
                            if (s3.substr(s3.find_first_not_of(" "), 7) == "<title>") {
                                newInData.setDetails(
                                        newInData.getDetails() + s3.substr(s3.find_first_not_of(" ") + 7,
                                                                           (s3.substr(
                                                                                   s3.find_first_not_of(
                                                                                           " "),
                                                                                   s3.find_last_not_of(
                                                                                           " "))).length()
                                                                           - 15) + " ");
                                //Break to move on to next description
                                break;
                            }
                        }
                    } else if (s2.substr(s2.find_first_not_of(" "), 9) == "<Hit_num>") {
                        //The extraction of data from nucleic acid BLAST results
                        for (int k = j + 1; k < sVec.size(); k++) {
                            std::string s3 = sVec.at(k);
                            //Check <id> tag to get the identifier
                            if (s3.substr(s3.find_first_not_of(" "), 8) == "<Hit_id>") {
                                //To avoid taking tags
                                std::string singleGeneDesc = s3.substr(s3.find_first_not_of(" ") + 8,
                                                                       (s3.substr(s3.find_first_not_of(" "),
                                                                                  s3.find_last_not_of(
                                                                                          " "))).length() -
                                                                       17);
                                //Add the description to the GeneData object
                                singleGeneDesc = singleGeneDesc.substr(0, singleGeneDesc.find_last_of("|"));
                                newInData.setDetails(newInData.getDetails() +
                                                     singleGeneDesc.substr(singleGeneDesc.find_last_of("|") + 1) +
                                                     " ");
                            }

                            //Check <title> tag for the additional identifier. (The scientific name etc.)
                            if (s3.substr(s3.find_first_not_of(" "), 9) == "<Hit_def>") {
                                newInData.setDetails(
                                        newInData.getDetails() + s3.substr(s3.find_first_not_of(" ") + 9,
                                                                           (s3.substr(
                                                                                   s3.find_first_not_of(
                                                                                           " "),
                                                                                   s3.find_last_not_of(
                                                                                           " "))).length()
                                                                           - 19) + " ");
                                //Break to move on to next description
                                break;
                            }
                        }
                    }
                }
                if (SupMorfi::getIdentifier(newInData.getDetails()).length() > 0 && newInData.getGene().length() > 0) {
                    //Make adding the new GeneData object to the output vector thread safe.
#pragma omp critical
                    {
                        inputData.push_back(newInData);
                    }
                }
            }
        }
        //Return the transformed data
        if (inputData.size() > 0) {
            return inputData;
        } else {
            throw IncorrectInpDataException();
        }
    } else {
        //Throw exception if an unknown file type
        throw UnknownInpFileException();
    }
}

// Get custom BLAST data in a custom file (not NCBI)
std::vector<GeneData> File2Obj::getInCustomBLASTData(std::string path, int thrdCnt) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

    //txt extraction
    if (ext == "txt") {

        //Read the entire file and convert it into a vector to support parallelism
        std::vector<std::string> sVec = getInpFileAsStrings(path);

//OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
        for (int i = 0; i < sVec.size(); i++) {
            //Check FASTA description start
            std::string s1 = sVec.at(i);
            if (s1.length() > 0 && s1[0] == '>') {
                //New GeneData Object
                GeneData newInData("", "");
                newInData.setDetails(s1.substr(1));

                int stopJ;
                //Iterate and add to description till 'Length' field or blank line
                for (int j = i; j < sVec.size(); j++) {
                    std::string s2 = sVec.at(j);
                    if ((s2.length() == 0) ||
                        (s2.length() > 0 && s2[0] == 'L' && s2[1] == 'e' && s2[2] == 'n')) {
                        stopJ = j;
                        break;
                    } else {
                        newInData.setDetails(newInData.getDetails() + s2);
                    }
                }

                //Temporarily hold the gene sequence
                std::string tempGene = "";
                //Temporarily hold the possible gene sequence part
                std::string tempLine = "";
                //Iterate through the FASTA sequences
                //Continue from after the gene's description
                for (int k = stopJ + 1; k < sVec.size(); k++) {
                    std::string s3 = sVec.at(k);
                    if (s3.length() > 0) {
                        //Stop if new description or file end component
                        if ((s3[0] == '>') || (s3[0] == ' ' && s3[1] == ' ' && s3[2] == 'D' && s3[3] == 'a' &&
                                               s3[4] == 't')) {
                            break;
                        }

                        //Separate the string line by spaces
                        std::istringstream buf(s3);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Take the latest string line longer than 10 characters
                        for (int q = 0; q < tokens.size(); q++) {
                            if (SupMorfi::isAlpha(tokens[q]) && tokens[q].length() >= 10) {
                                tempLine = tokens[q];
                            }
                        }

                    } else {
                        //When a blank space is found add the tempLine value to the tempGene and reset the tempLine value.
                        tempGene += tempLine;
                        tempLine = "";
                    }
                }
                //Set the resultant gene sequence
                newInData.setGene(tempGene);
                if (newInData.getGene().length() > 0) {
                    //Make adding the new GeneData object to the output vector thread safe.
#pragma omp critical
                    {
                        inputData.push_back(newInData);
                    }
                }
            }
        }
        //Return the transformed data
        if (inputData.size() > 0) {
            return inputData;
        } else {
            //If no matching GeneData(ble) components are found.
            throw IncorrectInpDataException();
        }

    } else {
        //Throw exception if an unknown file type
        throw UnknownInpFileException();
    }
}

// Get custom MSA(Multiple Sequence Alignment) data in a custom file (not NCBI)
std::vector<GeneData> File2Obj::getInCustomMSAData(std::string path) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

    //txt extraction
    if (ext == "txt") {
        //Input file stream
        std::ifstream infile;

        //Open file
        infile.open(path);

        //If there was no error in opening the file
        if (infile.good()) {
            int i = 0;
            //Boolean variable to inform if the first MSA segment where the GeneData objects are to be created are made
            //or not
            bool beginGene = true;

            //Temporary string to read the file
            std::string s;

            while (getline(infile, s)) {
                if (s != "") {
                    //If still in the first segment
                    if (beginGene) {
                        //Create new GeneData object to store the data
                        GeneData newInData("", "");

                        //Separate the string by spaces(" ")
                        std::istringstream buf(s);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Take the initial word (the NCBI or another ID to identify the gene) as the description
                        newInData.setDetails(tokens[0]);

                        //Add the later segments to the gene sequence
                        //*Spaces within the gene itself can be represented by spaces or other non-alphabetic characters
                        for (int l = 1; l < tokens.size(); l++) {
                            for (char c:tokens[l]) {
                                int cha = (int) c;
                                //If the character is alphabetic then add the charater to the gene
                                if (((cha <= (int) 'z' && cha >= (int) 'a') ||
                                     (cha <= (int) 'Z' && cha >= (int) 'A'))) {
                                    newInData.setGene(newInData.getGene() + c);
                                }
                            }
                        }
                        //Add the extracted gene to the output vector
                        inputData.push_back(newInData);
                    } else {
                        //Get the i th gene in the MSA segment
                        GeneData newInData = inputData.at(i);

                        //Take the initial word (the NCBI or another ID to identify the gene) as the description
                        std::istringstream buf(s);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Add the later segments to the gene sequence
                        //*Spaces within the gene itself can be represented by spaces or other non-alphabetic characters
                        for (int l = 1; l < tokens.size(); l++) {
                            for (char c:tokens[l]) {
                                int cha = (int) c;
                                if (((cha <= (int) 'z' && cha >= (int) 'a') ||
                                     (cha <= (int) 'Z' && cha >= (int) 'A'))) {
                                    newInData.setGene(newInData.getGene() + c);
                                }
                            }
                        }
                        //Update the gene
                        inputData.at(i) = newInData;
                    }
                    //Increment the number of alterations done for gene data
                    i += 1;
                } else {
                    //In the situation of a space
                    if (i != 0) {
                        //If genes have been altered and a space has been found set beginGene to false to indicate that
                        //the GeneData objects have been made
                        beginGene = false;
                    }
                    //Reset alterations. Used to identify the gene being altered in the stages other than the initial.
                    i = 0;
                }
            }
            //Return the transformed data
            if (inputData.size() > 0) {
                return inputData;
            } else {
                //Throw an exception if there are no matching data
                throw IncorrectInpDataException();
            }
        } else {
            //Throw error in file open
            throw FileInpException();
        }
    } else {
        //Throw exception if an unknown file type
        throw UnknownInpFileException();
    }
}

//The FASTA data file conversion function.
std::vector<GeneData> File2Obj::getInFASTAData(std::string path, int thrdCnt) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //Read the entire file and convert it into a vector to support parallelism
    std::vector<std::string> sVec = getInpFileAsStrings(path);

    //Read entire file, ( FASTA section by FASTA section )
    //OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
    for (int i = 0; i < sVec.size(); i++) {
        //Check FASTA description start
        std::string s1 = sVec.at(i);
        if (s1.length() > 0 && s1[0] == '>') {
            //New GeneData Object
            GeneData newInData("", "");
            //Add the first line to the description
            newInData.setDetails(s1.substr(1));

            //Iterate through the FASTA sequences
            for (int k = i + 1; k < sVec.size(); k++) {
                std::string s2 = sVec.at(k);
                if (s2.length() > 0) {
                    //Break the iteration if a new description has been found.
                    if ((s2[0] == '>')) {
                        break;
                    } else {
                        //Else take that as the input gene
                        newInData.setGene(newInData.getGene() + s2);
                    }
                }
            }
            //If a gene component has been found and is valid add it to the output vector
            if (newInData.getGene().length() > 0) {
                //Make adding the new GeneData object to the output vector thread safe.
#pragma omp critical
                {
                    inputData.push_back(newInData);
                }
            }
        }
    }
    //Return the transformed data
    if (inputData.size() > 0) {
        return inputData;
    } else {
        //If no matching GeneData(ble) components are found.
        throw IncorrectInpDataException();
    }
}

//The GCG data file conversion function.
std::vector<GeneData> File2Obj::getInGCGData(std::string path, int thrdCnt) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //Read the entire file and convert it into a vector to support parallelism
    std::vector<std::string> sVec = getInpFileAsStrings(path);

    //Read entire file, ( Gene by gene )
    //OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
    for (int i = 0; i < sVec.size(); i++) {
        //Check gene description start
        std::string s1 = sVec.at(i);
        if (s1.length() > 0 && s1.substr(0, 2) == "!!") {
            //New GeneData Object
            GeneData newInData("", "");

            //sD1 and sD2 would hold the descriptions of the gene.
            std::string sD1;
            std::string sD2;
            int i1, i2;

            //Take the first non empty string line
            for (i1 = i + 1; i1 < sVec.size(); i1++) {
                sD1 = sVec.at(i1);
                if (sD1.length() > 0) {
                    break;
                }
            }

            //Take the identifier from the second non empty sting line
            for (i2 = i1 + 1; i2 < sVec.size(); i2++) {
                std::string sD2Temp = sVec.at(i2);
                sD2 = sD2Temp.substr(0, sD2Temp.find_first_of(" "));
                if (sD2.length() > 0) {
                    break;
                }
            }
            //Make the gene description with the identifier and then the other
            newInData.setDetails(sD2 + " " + sD1);

            //Iterate through the gene sequences
            for (int k = i2 + 1; k < sVec.size(); k++) {
                std::string s2 = sVec.at(k);
                if (s2.length() > 0) {
                    //Break if a new gene is found
                    if ((s2.substr(0, 2) == "!!")) {
                        break;
                    } else {
                        //Split the line by spaces
                        std::istringstream buf(s2);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Add the space separated components except the first
                        for (int tokI = 1; tokI < tokens.size(); tokI++) {
                            newInData.setGene(newInData.getGene() + tokens[tokI]);
                        }
                    }
                }
            }
            //If a gene component has been found and is valid add it to the output vector
            if (newInData.getGene().length() > 0) {
                //Make adding the new GeneData object to the output vector thread safe.
#pragma omp critical
                {
                    inputData.push_back(newInData);
                }
            }
        }
    }
    //Return the transformed data
    if (inputData.size() > 0) {
        return inputData;
    } else {
        //If no matching GeneData(ble) components are found.
        throw IncorrectInpDataException();
    }
}

//The EMBL data file conversion function.
std::vector<GeneData> File2Obj::getInEMBLData(std::string path, int thrdCnt) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //Read the entire file and convert it into a vector to support parallelism
    std::vector<std::string> sVec = getInpFileAsStrings(path);

    //Read entire file, ( Gene by gene )
    //OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
    for (int i = 0; i < sVec.size(); i++) {
        //Check description start
        std::string s1 = sVec.at(i);
        if (s1.length() > 0 && s1.substr(0, 2) == "ID") {
            //New GeneData Object
            GeneData newInData("", "");
            //Get the gene ID
            std::string temID = s1.substr(3);
            newInData.setDetails(temID.substr(temID.find_first_not_of(" "),
                                              temID.find_first_of(';') - temID.find_first_not_of(" ")));

            //Find and locate the description components other than the ID
            std::string sD1;
            int i1;
            for (i1 = i + 1; i1 < sVec.size(); i1++) {
                sD1 = sVec.at(i1);
                if (sD1.substr(0, 2) == "DE") {
                    std::string sD2 = sD1.substr(3);
                    newInData.setDetails(newInData.getDetails() + " " + sD2.substr(sD2.find_first_not_of(" ")));
                } else if (sD1.length() > 0 && sD1.substr(0, 2) == "SQ") {
                    //Break if the sequence component arises
                    break;
                }
            }

            //Iterate through the gene sequences
            for (int k = i1 + 1; k < sVec.size(); k++) {
                std::string s2 = sVec.at(k);
                if (s2.length() > 0) {
                    //Break if a new gene is found
                    if ((s2.substr(0, 2) == "//")) {
                        break;
                    } else {
                        //Split the line by spaces
                        std::istringstream buf(s2);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Add the space separated components except the first
                        for (int tokI = 0; tokI < tokens.size() - 1; tokI++) {
                            newInData.setGene(newInData.getGene() + tokens[tokI]);
                        }
                    }
                }
            }
            //If a gene component has been found and is valid add it to the output vector
            if (newInData.getGene().length() > 0) {
                //Make adding the new GeneData object to the output vector thread safe.
#pragma omp critical
                {
                    inputData.push_back(newInData);
                }
            }
        }
    }
    //Return the transformed data
    if (inputData.size() > 0) {
        return inputData;
    } else {
        //If no matching GeneData(ble) components are found.
        throw IncorrectInpDataException();
    }
}

//The Genbank data file conversion function.
std::vector<GeneData> File2Obj::getInGenebankData(std::string path, int thrdCnt) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //Read the entire file and convert it into a vector to support parallelism
    std::vector<std::string> sVec = getInpFileAsStrings(path);

    //Read entire file, ( Gene by gene )
    //OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
    for (int i = 0; i < sVec.size(); i++) {
        //Check FASTA description start
        std::string s1 = sVec.at(i);
        if (s1.length() > 0 && s1.substr(s1.find_first_not_of(" "), 5) == "LOCUS") {
            //New GeneData Object
            GeneData newInData("", "");

            //Get the gene ID
            std::string temID1 = s1.substr(s1.find_first_not_of(" ") + 5);
            std::string temID2 = temID1.substr(temID1.find_first_not_of(" "));
            newInData.setDetails(temID2.substr(0, temID2.find_first_of(" ")));

            //Find and locate the description components other than the ID
            std::string sD1;
            int i1;
            for (i1 = i + 1; i1 < sVec.size(); i1++) {
                sD1 = sVec.at(i1);
                if (sD1.substr(sD1.find_first_not_of(" "), 10) == "DEFINITION") {
                    std::string temDes1 = sD1.substr(sD1.find_first_not_of(" ") + 10);

                    newInData.setDetails(
                            newInData.getDetails() + " " + temID1.substr(temID1.find_first_not_of(" ")));
                } else if (sD1.length() > 0 && sD1.substr(sD1.find_first_not_of(" "), 6) == "ORIGIN") {
                    //Break if the sequence component arises
                    break;
                }
            }

            //Iterate through the gene sequences
            for (int k = i1 + 1; k < sVec.size(); k++) {
                std::string s2 = sVec.at(k);
                if (s2.length() > 0) {
                    //Break if a new gene is found
                    if ((s2.substr(s2.find_first_not_of(" "), 2) == "//")) {
                        break;
                    } else {
                        //Split the line by spaces
                        std::istringstream buf(s2);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Add the space separated components except the first
                        for (int tokI = 1; tokI < tokens.size(); tokI++) {
                            newInData.setGene(newInData.getGene() + tokens[tokI]);
                        }
                    }
                }
            }
            //If a gene component has been found and is valid add it to the output vector
            if (newInData.getGene().length() > 0) {
                //Make adding the new GeneData object to the output vector thread safe.
#pragma omp critical
                {
                    inputData.push_back(newInData);
                }
            }
        }
    }
    //Return the transformed data
    if (inputData.size() > 0) {
        return inputData;
    } else {
        //If no matching GeneData(ble) components are found.
        throw IncorrectInpDataException();
    }
}

//The PIR data file conversion function.
std::vector<GeneData> File2Obj::getInPIRData(std::string path, int thrdCnt) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //Read the entire file and convert it into a vector to support parallelism
    std::vector<std::string> sVec = getInpFileAsStrings(path);

    //Read entire file, ( Gene by gene )
    //OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
    for (int i = 0; i < sVec.size(); i++) {
        //Check FASTA description start
        std::string s1 = sVec.at(i);
        if (s1.length() > 1 && s1.substr(s1.find_first_not_of(" "))[0] == '>') {

            //New GeneData Object
            GeneData newInData("", "");
            //Get the gene ID
            newInData.setDetails(s1.substr(s1.find_first_of(";") + 1));

            //Find and locate the description components other than the ID
            std::string sD1 = sVec.at(i + 1);
            std::string sD2 = sD1.substr(0, sD1.find_last_of(","));
            newInData.setDetails(newInData.getDetails() + " " + sD2.substr(0, sD2.find_last_not_of(" ")));

            //Iterate through the gene sequences
            for (int k = i + 2; k < sVec.size(); k++) {
                std::string s2 = sVec.at(k);
                if (s2.length() > 0) {
                    //Break if a new gene is found
                    if ((s2.substr(s2.find_first_not_of(" "), 1) == ">" || s2.length() == 0)) {
                        break;
                    } else {
                        //Split the line by spaces
                        std::istringstream buf(s2);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Add the space separated components
                        for (int tokI = 0; tokI < tokens.size(); tokI++) {
                            newInData.setGene(newInData.getGene() + SupMorfi::extractAlpha(tokens[tokI]));
                        }
                    }
                }
            }
            //If a gene component has been found and is valid add it to the output vector
            if (newInData.getGene().length() > 0) {
                //Make adding the new GeneData object to the output vector thread safe.
#pragma omp critical
                {
                    inputData.push_back(newInData);
                }
            }
        }
    }
    //Return the transformed data
    if (inputData.size() > 0) {
        return inputData;
    } else {
        //If no matching GeneData(ble) components are found.
        throw IncorrectInpDataException();
    }
}

//The PHYLI data file conversion function.
std::vector<GeneData> File2Obj::getInPHYLIPData(std::string path, int thrdCnt) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //Read the entire file and convert it into a vector to support parallelism
    std::vector<std::string> sVec = getInpFileAsStrings(path);

    //Read entire file, ( Gene by gene )
    //OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
    for (int i = 1; i < sVec.size(); i++) {
        //Check FASTA description start
        std::string s1 = sVec.at(i);
        if (s1.length() > 0 && s1[0] != ' ') {

            //New GeneData Object
            GeneData newInData("", "");

            //Find the length of the ID by looking at the string line after the first
            std::string sD1 = sVec.at(i + 1);
            int spaceToGene = sD1.find_first_not_of(" ");

            //Get the gene ID
            std::string sD2 = s1.substr(0, spaceToGene);
            newInData.setDetails(sD2);

            //Data without the ID
            std::string sG = s1.substr(spaceToGene);

            //Split the remaining line by spaces
            std::istringstream bufTem(sG);
            std::istream_iterator<std::string> begTem(bufTem), endTem;
            std::vector<std::string> tokensTem(begTem, endTem);

            //Add the gene sequence in the first line
            for (int tokITem = 0; tokITem < tokensTem.size(); tokITem++) {
                newInData.setGene(newInData.getGene() + tokensTem[tokITem]);
            }

            //Iterate through the gene sequences
            for (int k = i + 1; k < sVec.size(); k++) {
                std::string s2 = sVec.at(k);
                if (s2.length() > 0) {
                    //Break if a new gene is found
                    if ((s2[0] != ' ' || s2.length() == 0)) {
                        break;
                    } else {
                        //Split the line by spaces
                        std::istringstream buf(s2);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Add the space separated components
                        for (int tokI = 0; tokI < tokens.size(); tokI++) {
                            newInData.setGene(newInData.getGene() + SupMorfi::extractAlpha(tokens[tokI]));
                        }
                    }
                }
            }
            //If a gene component has been found and is valid add it to the output vector
            if (newInData.getGene().length() > 0) {
                //Make adding the new GeneData object to the output vector thread safe.
#pragma omp critical
                {
                    inputData.push_back(newInData);
                }
            }
        }
    }
    //Return the transformed data
    if (inputData.size() > 0) {
        return inputData;
    } else {
        //If no matching GeneData(ble) components are found.
        throw IncorrectInpDataException();
    }
}

//The Clustral data file (A Multiple Sequence Alignment) conversion function.
std::vector<GeneData> File2Obj::getInClustralData(std::string path) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

    //txt extraction
    if (ext == "txt") {
        //Input file stream
        std::ifstream infile;

        //Open file
        infile.open(path);

        //If there was no error in opening the file
        if (infile.good()) {
            int i = 0;
            //Boolean variable to inform if the first MSA segment where the GeneData objects are to be created are made
            //or not
            bool beginGene = true;

            //Temporary string to read the file
            std::string s;

            //Skip the first line. (Mentions the file type)
            getline(infile, s);


            while (getline(infile, s)) {
                if (s.length() > 0 && s[0] != ' ') {
                    //If still in the first segment
                    if (beginGene) {
                        //Create new GeneData object to store the data
                        GeneData newInData("", "");

                        //Separate the string by spaces(" ")
                        std::istringstream buf(s);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //The identifier
                        std::string sD1 = tokens[0].substr(tokens[0].find_last_of('|') + 1);
                        //The gene description other than the identifier
                        std::string sD2 = tokens[0].substr(0, tokens[0].find_last_of('|'));

                        //Set the details of the gene
                        newInData.setDetails(sD1 + " " + sD2);

                        //Add the later segments to the gene sequence
                        //*Spaces within the gene itself can be represented by spaces or other non-alphabetic characters
                        for (int l = 1; l < tokens.size(); l++) {
                            for (char c:tokens[l]) {
                                int cha = (int) c;
                                if (((cha <= (int) 'z' && cha >= (int) 'a') ||
                                     (cha <= (int) 'Z' && cha >= (int) 'A'))) {
                                    newInData.setGene(newInData.getGene() + c);
                                }
                            }
                        }
                        inputData.push_back(newInData);
                    } else {
                        //Get the i th gene in the segment
                        GeneData newInData = inputData.at(i);

                        //Separate the string by spaces(" ")
                        std::istringstream buf(s);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Add the later segments to the gene sequence
                        //*Spaces within the gene itself can be represented by spaces or other non-alphabetic characters
                        for (int l = 1; l < tokens.size(); l++) {
                            for (char c:tokens[l]) {
                                int cha = (int) c;
                                if (((cha <= (int) 'z' && cha >= (int) 'a') ||
                                     (cha <= (int) 'Z' && cha >= (int) 'A'))) {
                                    newInData.setGene(newInData.getGene() + c);
                                }
                            }
                        }
                        //Update the gene
                        inputData.at(i) = newInData;
                    }
                    //Increment the number of alterations done for gene data
                    i += 1;
                } else {
                    //In the situation of a space
                    if (i != 0) {
                        //If genes have been altered and a space has been found set beginGene to false to indicate that
                        //the GeneData objects have been made
                        beginGene = false;
                    }
                    //Reset alterations. Used to identify the gene being altered in the stages other than the initial.
                    i = 0;
                }
            }
            //Return the transformed data
            if (inputData.size() > 0) {
                return inputData;
            } else {
                //Throw an exception if there are no matching data
                throw IncorrectInpDataException();
            }
        } else {
            //Throw error in file open
            throw FileInpException();
        }
    } else {
        //Throw exception if an unknown file type
        throw UnknownInpFileException();
    }
}

//The MSF data file conversion function.
std::vector<GeneData> File2Obj::getInMSFData(std::string path) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

    //txt extraction
    if (ext == "txt") {
        //Input file stream
        std::ifstream infile;

        //Open file
        infile.open(path);

        //If there was no error in opening the file
        if (infile.good()) {
            int i = 0;
            //Boolean variable to inform if the first MSA segment where the GeneData objects are to be created are made
            //or not
            bool beginGene = true;

            //Temporary string to read the file
            std::string s;

            //Skip till the gene data component. (Skip the summaries of the genes)
            while (getline(infile, s)) {
                if (s.length() > 0 && s.substr(0, 2) == "//") {
                    break;
                }
            }

            while (getline(infile, s)) {
                if (s.length() > 0 && s[0] != ' ') {
                    //If still in the first segment
                    if (beginGene) {
                        //Create new GeneData object to store the data
                        GeneData newInData("", "");

                        //Separate the string by spaces(" ")
                        std::istringstream buf(s);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        newInData.setDetails(tokens[0]);

                        //Add the later segments to the gene sequence
                        //*Spaces within the gene itself can be represented by spaces or other non-alphabetic characters
                        for (int l = 1; l < tokens.size(); l++) {
                            for (char c:tokens[l]) {
                                int cha = (int) c;
                                if (((cha <= (int) 'z' && cha >= (int) 'a') ||
                                     (cha <= (int) 'Z' && cha >= (int) 'A'))) {
                                    newInData.setGene(newInData.getGene() + c);
                                }
                            }
                        }
                        inputData.push_back(newInData);
                    } else {
                        //Get the i th gene in the segment
                        GeneData newInData = inputData.at(i);

                        //Separate the string by spaces(" ")
                        std::istringstream buf(s);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Add the later segments to the gene sequence
                        //*Spaces within the gene itself can be represented by spaces or other non-alphabetic characters
                        for (int l = 1; l < tokens.size(); l++) {
                            for (char c:tokens[l]) {
                                int cha = (int) c;
                                if (((cha <= (int) 'z' && cha >= (int) 'a') ||
                                     (cha <= (int) 'Z' && cha >= (int) 'A'))) {
                                    newInData.setGene(newInData.getGene() + c);
                                }
                            }
                        }
                        //Update the gene
                        inputData.at(i) = newInData;
                    }
                    //Increment the number of alterations done for gene data
                    i += 1;
                } else {
                    //In the situation of a space
                    if (i != 0) {
                        //If genes have been altered and a space has been found set beginGene to false to indicate that
                        //the GeneData objects have been made
                        beginGene = false;
                    }
                    //Reset alterations. Used to identify the gene being altered in the stages other than the initial.
                    i = 0;
                }
            }
            //Return the transformed data
            if (inputData.size() > 0) {
                return inputData;
            } else {
                //Throw an exception if there are no matching data
                throw IncorrectInpDataException();
            }
        } else {
            //Throw error in file open
            throw FileInpException();
        }
    } else {
        //Throw exception if an unknown file type
        throw UnknownInpFileException();
    }
}

//The FASTA Report data file conversion function.
std::vector<GeneData> File2Obj::getInFASTARepData(std::string path, int thrdCnt) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //Read the entire file and convert it into a vector to support parallelism
    std::vector<std::string> sVec = getInpFileAsStrings(path);

    //Read entire file, ( FASTA section by FASTA section )
    //OpenMP line to read the file content in "thrdCnt" number of threads.
#pragma omp parallel for num_threads(thrdCnt)
    for (int i = 1; i < sVec.size(); i++) {
        //Check FASTA description start
        std::string s1 = sVec.at(i);
        if (s1.length() > 0 && s1.substr(0, 2) == ">>") {

            //New GeneData Object
            GeneData newInData("", "");

            //Get the gene description
            newInData.setDetails(s1.substr(s1.find_first_of(" ")));
            //Get the identifier for the gene in the report
            std::string sLoc = s1.substr(2, s1.find_first_of(":") - 2);

            //Iterate through the gene sequences
            for (int k = i + 1; k < sVec.size(); k++) {
                std::string s2 = sVec.at(k);
                if (s2.length() > 0) {
                    //Break if a new gene is found
                    if (s2.substr(0, 2) == ">>") {
                        break;
                    } else {
                        //Split the line by spaces
                        std::istringstream buf(s2);
                        std::istream_iterator<std::string> beg(buf), end;
                        std::vector<std::string> tokens(beg, end);

                        //Add the gene component so long as the identifier matches
                        if (tokens[0] == sLoc) {
                            newInData.setGene(newInData.getGene() + tokens[1]);

                        }
                    }
                }
            }
            //If a gene component has been found and is valid add it to the output vector
            if (newInData.getGene().length() > 0) {
                //Make adding the new GeneData object to the output vector thread safe.
#pragma omp critical
                {
                    inputData.push_back(newInData);
                }
            }
        }
    }
    //Return the transformed data
    if (inputData.size() > 0) {
        return inputData;
    } else {
        //If no matching GeneData(ble) components are found.
        throw IncorrectInpDataException();
    }
}