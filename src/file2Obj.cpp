#include "file2Obj.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include <memory>

//#include <allheaders.h> // leptonica main header for image io
//#include <baseapi.h> // tesseract main header

#if defined(_OPENMP)

#include <omp.h>

#endif

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

//Incorrect input file type exception
class UnknownInpFileException : public MorfiFile2ObjException {
public:
    virtual const char *what() const throw() {
        return "Unsupported input file format";
    }
};



//Conversion functions

//Get the inputs form a NCBI BLAST result
std::vector<GeneData> File2Obj::getInNCBIBLASTData(std::string path) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

    //txt extraction
    if (ext == "txt") {

//Check if OpenMP is supported
// If OpenMP is supported
#if defined(_OPENMP)
        //Input file stream
        std::ifstream infile;

        //Open file
        infile.open(path);

        //If there was no error in opening the file
        if (infile.good()) {
            //String to hold every file-line read
            std::vector<std::string> sVec;

            //Temporary storage string
            std::string sTemp;

            // Add all the content in the file to sVec (To incorporate parallelism when processing the data)
            while (getline(infile, sTemp)) {
                sVec.push_back(sTemp);
            }

            //Close and clear input file stream once done
            infile.close();
            infile.clear();

            //Read entire file, ( FASTA section by FASTA section )
#pragma omp parallel num_threads(10)
            {
#pragma omp for
                for (int i = 0; i < sVec.size(); i++) {
                    //Check FASTA description start
                    std::string s1 = sVec.at(i);
                    if (s1.length() > 0 && s1[0] == '>') {
                        int pointerFASTA = 1;
                        //New GeneData Object
                        GeneData newInData("", "");
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
                                    break;
                                }
                                pointerFASTA = currentPointer;

                            } else if (s3.length() > 0 && ((s3[0] == '>') ||
                                                           (s3[0] == ' ' && s3[1] == ' ' && s3[2] == 'D' &&
                                                            s3[3] == 'a' &&
                                                            s3[4] == 't'))) {
                                //Add the GeneData to the inputData set once a new FASTA description is discovered and recognized
                                if (newInData.getGene().length() > 0) {
                                    inputData.push_back(newInData);
                                }
                                break;
                            }
                        }
                    }
                }
            }
//If OpenMP is not supported
#else
            //Input file stream
            std::ifstream infile;

            //Open file
            infile.open(path);

            //If there was no error in opening the file
            if (infile.good()) {
                //String to hold every file-line read
                std::string s;
                //Current file pointer/position
                std::streampos pos = 0;

                //Read entire file, ( FASTA section by FASTA section )
                while (getline(infile, s)) {
                    //Check FASTA description start
                    if (s.length() > 0 && s[0] == '>') {
                        int pointerFASTA = 1;
                        //New GeneData Object
                        GeneData newInData("", "");
                        newInData.setDetails(s.substr(1));

                        //Iterate and add to description till 'Length' field
                        while (getline(infile, s)) {
                            if (s.length() > 0 && s[0] == 'L' && s[1] == 'e' && s[2] == 'n') {
                                break;
                            } else {
                                newInData.setDetails(newInData.getDetails() + s);
                            }
                        }

                        //Iterate through the FASTA sequences
                        while (getline(infile, s)) {
                            //Check if subject FASTA
                            if (s.length() > 0 && s[0] == 'S' && s[1] == 'b' && s[2] == 'j' && s[3] == 'c' && s[4] == 't') {

                                //String stream to iterate through a line and find the necessary sections
                                std::istringstream buf(s);
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
                                    break;
                                }
                                pointerFASTA = currentPointer;

                            } else if (s.length() > 0 && ((s[0] == '>') ||
                                                          (s[0] == ' ' && s[1] == ' ' && s[2] == 'D' && s[3] == 'a' &&
                                                           s[4] == 't'))) {
                                //Add the GeneData to the inputData set once a new FASTA description is discovered and recognized
                                if (newInData.getGene().length() > 0) {
                                    inputData.push_back(newInData);
                                }

                                //Reposition to the previous line
                                infile.clear();
                                infile.seekg(pos);
                                break;
                            }
                            //Get the position of the current(already read) line and keep as to reposition as necessary
                            pos = infile.tellg();
                        }
                    }
                }

                //Close and clear input file stream once done
                infile.close();
                infile.clear();
#endif

            //Return the transformed data
            return inputData;
        } else {
            //Throw error in file open
            throw FileInpException();
        }
        //json extraction
    } else if (ext == "json") {
        //Variable to hold all json values
        Json::Value root;

        //Get file and parse the values in it as json and assign to the holder variable
        std::ifstream infile(path);
        infile >> root;

        //Array that the BLAST results iterate over
        const Json::Value &res = root["BlastOutput2"]["report"]["results"]["search"]["hits"];

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
                geneDisc += (exrt.substr(exrt.find_last_of('|') + 1) + " ");
                //Add the title to the gene sequence description
                geneDisc += (gene["title"].asString() + " ");
            }
            //Remove tailing whitespace
            newInData.setDetails(geneDisc.substr(0, geneDisc.length() - 1));
            //Set the FASTA sequence
            newInData.setGene(geneSet["hsps"][0]["hseq"].asString());
            //Add the gene data to the vector set
            inputData.push_back(newInData);
        }

        //Return full data set
        return inputData;
    } else if (ext == "xml") {
#if defined(_OPENMP)
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
            throw FileInpException();
        }

        //Read entire file, ( FASTA section by FASTA section )
#pragma omp parallel num_threads(10)
        {
#pragma omp for
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
                        if (s2.substr(s2.find_first_not_of(" "), 6) == "<hseq>") {
                            //Get substring to avoid taking xml tags
                            newInData.setGene(
                                    s2.substr(s2.find_first_not_of(" ") + 6, s2.substr(s2.find_first_not_of(" ")).
                                            length() - 13));
                            break;
                        }

                        //A description of a gene sequence, iterate within to get the NCBI identifier and the additional
                        //description. AVOID LEADING WHITE SPACES / REMEMBER TO TRIM THEM!!
                        if (s2.substr(s2.find_first_not_of(" "), s2.find_last_of(" ")) == "<HitDescr>") {
                            for (int k = j + 1; k < sVec.size(); k++) {
                                std::string s3 = sVec.at(k);
                                //Check <id> tag to get the identifier
                                if (s3.substr(s3.find_first_not_of(" "), 4) == "<id>") {
                                    //To avoid taking tags,
                                    std::string singleGeneDesc = s3.substr(s3.find_first_not_of(" ") + 4,
                                                                           (s3.substr(s3.find_first_not_of(" "),
                                                                                      s3.find_last_not_of(
                                                                                              " "))).length() -
                                                                           9);
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
                        }
                    }
                    //Remove tailing white space.
                    //newInData.setDetails(newInData.getDetails().substr(0,newInData.getDetails().length()-1));
                    //Add gathered description data to the final vector set
                    inputData.push_back(newInData);
                }
            }
        }
#else
        std::ifstream infile;

        //Open file
        infile.open(path);

        //If there was no error in opening the file
        if (infile.good()) {
            //String to hold every file-line read
            std::string s;
            //Current file pointer/position
            std::streampos pos = 0;

            //Read entire file, ( FASTA section by FASTA section )
            while (getline(infile, s)) {

                //Start of a gene sequence used for BLAST
                if (s.substr(s.find_first_not_of(" "), s.find_last_of(" ")) == "<Hit>") {
                    //New GeneData object
                    GeneData newInData("", "");
                    //To hold the description of the gene sequence being processed
                    std::string genDisc = "";

                    //iterate through all the gene sequence descriptions and FASTA sequence
                    while (getline(infile, s)) {
                        //If the FASTA sequence is located then add to GeneData and break from loop
                        if (s.substr(s.find_first_not_of(" "), 6) == "<hseq>") {
                            //Get substing to avoid taking xml tags
                            newInData.setGene(s.substr(s.find_first_not_of(" ") + 6, s.substr(s.find_first_not_of(" ")).
                                    length() - 13));
                            break;
                        }

                        //A description of a gene sequence, iterate within to get the NCBI identifier and the additional
                        //description. AVOID LEADING WHITE SPACES / REMEMBER TO TRIM THEM!!
                        if (s.substr(s.find_first_not_of(" "), s.find_last_of(" ")) == "<HitDescr>") {
                            while (getline(infile, s)) {
                                //Check <id> tag to get the identifier
                                if (s.substr(s.find_first_not_of(" "), 4) == "<id>") {
                                    //To avoid taking tags,
                                    std::string singleGeneDesc = s.substr(s.find_first_not_of(" ") + 4,
                                                                          (s.substr(s.find_first_not_of(" "),
                                                                                    s.find_last_not_of(" "))).length() -
                                                                          9);
                                    singleGeneDesc = singleGeneDesc.substr(0, singleGeneDesc.find_last_of("|"));
                                    newInData.setDetails(newInData.getDetails() +
                                                         singleGeneDesc.substr(singleGeneDesc.find_last_of("|") + 1) +
                                                         " ");
                                }

                                //Check <title> tag for the additional identifier. (The scientific name etc.)
                                if (s.substr(s.find_first_not_of(" "), 7) == "<title>") {
                                    newInData.setDetails(newInData.getDetails() + s.substr(s.find_first_not_of(" ") + 7,
                                                                                           (s.substr(
                                                                                                   s.find_first_not_of(
                                                                                                           " "),
                                                                                                   s.find_last_not_of(
                                                                                                           " "))).length()
                                                                                           - 15) + " ");
                                    //Break to move on to next description
                                    break;
                                }
                            }
                        }
                    }
                    //Remove tailing white space.
                    //newInData.setDetails(newInData.getDetails().substr(0,newInData.getDetails().length()-1));
                    //Add gathered description data to the final vector set
                    inputData.push_back(newInData);
                }
            }
        } else {
            throw FileInpException();
        }
#endif
        return inputData;
    } else {
        //Throw exception if an unknown file type
        throw UnknownInpFileException();
    }
}

//Get the inputs form a NCBI BLAST result
std::vector<GeneData> File2Obj::getInCustomBLASTData(std::string path) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

    //txt extraction
    if (ext == "txt") {

//Check if OpenMP is supported
// If OpenMP is supported
#if defined(_OPENMP)
        //Input file stream
        std::ifstream infile;

        //Open file
        infile.open(path);

        //If there was no error in opening the file
        if (infile.good()) {
            //String to hold every file-line read
            std::vector<std::string> sVec;

            //Temporary storage string
            std::string sTemp;

            // Add all the content in the file to sVec (To incorporate parallelism when processing the data)
            while (getline(infile, sTemp)) {
                sVec.push_back(sTemp);
            }

            //Close and clear input file stream once done
            infile.close();
            infile.clear();

            //Read entire file, ( FASTA section by FASTA section )
#pragma omp parallel num_threads(10)
            {
#pragma omp for
                for (int i = 0; i < sVec.size(); i++) {
                    //Check FASTA description start
                    std::string s1 = sVec.at(i);
                    if (s1.length() > 0 && s1[0] == '>') {
                        //New GeneData Object
                        GeneData newInData("", "");
                        newInData.setDetails(s1.substr(1));

                        int stopJ;
                        //Iterate and add to description till 'Length' field
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

                        std::string tempGene = "";
                        std::string tempLine = "";
                        //Iterate through the FASTA sequences
                        for (int k = stopJ + 1; k < sVec.size(); k++) {
                            std::string s3 = sVec.at(k);
                            if (s3.length() > 0) {
                                if ((s3[0] == '>') || (s3[0] == ' ' && s3[1] == ' ' && s3[2] == 'D' && s3[3] == 'a' &&
                                                       s3[4] == 't')) {
                                    newInData.setGene(tempGene);
                                    break;
                                }

                                std::istringstream buf(s3);
                                std::istream_iterator<std::string> beg(buf), end;

                                std::vector<std::string> tokens(beg, end);

                                for (int l = 0; l < tokens.size(); l++) {
                                    if (SupMorfi::isAlpha(tokens[l]) && tokens[l].length() >= 30) {
                                        tempLine = tokens[l];
                                    }
                                }

                            } else {
                                tempGene += tempLine;
                            }
                        }
                        if (newInData.getGene().length() > 0) {
                            inputData.push_back(newInData);
                        }
                    }
                }
            }
//If OpenMP is not supported
#else
            //Input file stream
            std::ifstream infile;

            //Open file
            infile.open(path);

            //If there was no error in opening the file
            if (infile.good()) {
                //String to hold every file-line read
                std::string s;
                //Current file pointer/position
                std::streampos pos = 0;

                //Read entire file, ( FASTA section by FASTA section )
                while (getline(infile, s)) {
                    //Check FASTA description start
                    if (s.length() > 0 && s[0] == '>') {
                        int pointerFASTA = 1;
                        //New GeneData Object
                        GeneData newInData("", "");
                        newInData.setDetails(s.substr(1));

                        //Iterate and add to description till 'Length' field
                        while (getline(infile, s)) {
                            if ((s.length() == 0) || (s.length() > 0 && s[0] == 'L' && s[1] == 'e' && s[2] == 'n'))  {
                                break;
                            } else {
                                newInData.setDetails(newInData.getDetails() + s);
                            }
                        }

                        std::string tempGene="";
                        std::string tempLine="";
                        //Iterate through the FASTA sequences
                        while (getline(infile, s)) {
                            //Check if subject FASTA
                            if (s.length() > 0) {
                                if ((s[0] == '>') || (s[0] == ' ' && s[1] == ' ' && s[2] == 'D' && s[3] == 'a' &&
                                     s[4] == 't')){
                                    newInData.setGene(tempGene);
                                    //Add the GeneData to the inputData set once a new FASTA description is discovered and recognized
                                    if (newInData.getGene().length() > 0) {
                                        inputData.push_back(newInData);
                                    }

                                    //Reposition to the previous line
                                    infile.clear();
                                    infile.seekg(pos);
                                    break;
                                }
                                //String stream to iterate through a line and find the necessary sections
                                std::istringstream buf(s);
                                std::istream_iterator<std::string> beg(buf), end;

                                std::vector<std::string> tokens(beg, end);

                                //Add to the FASTA sequence
                                //Check pointer in FASTA, in case a description is missing.

                                for (int l = 0; l < tokens.size(); l++) {
                                    if (SupFile2ObjFnc::isAlpha(tokens[l]) && tokens[l].length()>=30){
                                        tempLine=tokens[l];
                                    }
                                }
                            } else {
                                tempGene+=tempLine;
                            }
                            //Get the position of the current(already read) line and keep as to reposition as necessary
                            pos = infile.tellg();
                        }
                    }
                }

                //Close and clear input file stream once done
                infile.close();
                infile.clear();
#endif

            //Return the transformed data
            return inputData;
        } else {
            //Throw error in file open
            throw FileInpException();
        }
    } else {
        //Throw exception if an unknown file type
        throw UnknownInpFileException();
    }
}

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
            bool beginGene = true;
            std::string s;

            while (getline(infile, s)) {
//                std::cout << s <<std::endl;
                if (s != "") {
                    if (beginGene) {
                        GeneData newInData("", "");

                        std::istringstream buf(s);
                        std::istream_iterator<std::string> beg(buf), end;

                        std::vector<std::string> tokens(beg, end);

                        newInData.setDetails(tokens[0]);

                        for (int l = 1; l < tokens.size(); l++) {
                            for (char c:tokens[l]) {
                                int cha = (int) c;
                                if (((cha <= (int) 'z' && cha >= (int) 'a') ||
                                     (cha <= (int) 'Z' && cha >= (int) 'A'))) {
                                    newInData.setGene(newInData.getGene() + c);
                                }
                            }
                        }
//                        std::cout<< newInData.getGene()<<std::endl;
                        inputData.push_back(newInData);
                    } else {
                        GeneData newInData = inputData.at(i);
//                        std::cout<< newInData.getDetails()<<std::endl;

                        std::istringstream buf(s);
                        std::istream_iterator<std::string> beg(buf), end;

                        std::vector<std::string> tokens(beg, end);

                        for (int l = 1; l < tokens.size(); l++) {
                            for (char c:tokens[l]) {
                                int cha = (int) c;
                                if (((cha <= (int) 'z' && cha >= (int) 'a') ||
                                     (cha <= (int) 'Z' && cha >= (int) 'A'))) {
                                    newInData.setGene(newInData.getGene() + c);
                                }
                            }
                        }
                        inputData.at(i) = newInData;
                    }
                    i += 1;
                } else {
                    if (i != 0) {
                        beginGene = false;
                    }
                    i = 0;
                }
            }
            return inputData;
        } else {
            //Throw error in file open
            throw FileInpException();
        }
    } else {
//Throw exception if an unknown file type
        throw UnknownInpFileException();
    }
}


//Get the inputs form a NCBI BLAST result
std::vector<GeneData> File2Obj::getInFASTAData(std::string path) {
    //Vector set to hold converted input data set
    std::vector<GeneData> inputData;

    //Get the extension type to perform corresponding extraction
    std::string ext = SupMorfi::getExtension(path);

//Check if OpenMP is supported
// If OpenMP is supported
#if defined(_OPENMP)
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //If there was no error in opening the file
    if (infile.good()) {
        //String to hold every file-line read
        std::vector<std::string> sVec;

        //Temporary storage string
        std::string sTemp;

        // Add all the content in the file to sVec (To incorporate parallelism when processing the data)
        while (getline(infile, sTemp)) {
            sVec.push_back(sTemp);
        }

        //Close and clear input file stream once done
        infile.close();
        infile.clear();

        //Read entire file, ( FASTA section by FASTA section )
#pragma omp parallel num_threads(10)
        {
#pragma omp for
            for (int i = 0; i < sVec.size(); i++) {
                //Check FASTA description start
                std::string s1 = sVec.at(i);
                if (s1.length() > 0 && s1[0] == '>') {
                    //New GeneData Object
                    GeneData newInData("", "");
                    newInData.setDetails(s1.substr(1));

                    //Iterate through the FASTA sequences
                    for (int k = i + 1; k < sVec.size(); k++) {
                        std::string s2 = sVec.at(k);
                        if (s2.length() > 0) {
                            if ((s2[0] == '>')) {
                                break;
                            } else {
                                newInData.setGene(newInData.getGene() + s2);
                            }
                        }
                    }
                    if (newInData.getGene().length() > 0) {
                        inputData.push_back(newInData);
                    }
                }
            }
        }
//If OpenMP is not supported
#else
        //Input file stream
        std::ifstream infile;

        //Open file
        infile.open(path);

        //If there was no error in opening the file
        if (infile.good()) {
            //String to hold every file-line read
            std::string s;
            //Current file pointer/position
            std::streampos pos = 0;

            //Read entire file, ( FASTA section by FASTA section )
            while (getline(infile, s)) {
                //Check FASTA description start
                if (s.length() > 0 && s[0] == '>') {
                    int pointerFASTA = 1;
                    //New GeneData Object
                    GeneData newInData("", "");
                    newInData.setDetails(s.substr(1));

                    //Iterate through the FASTA sequences
                    while (getline(infile, s)) {
                        //Check if subject FASTA
                        if (s.length() > 0) {
                            if ((s[0] == '>')){
                                //Reposition to the previous line
                                infile.clear();
                                infile.seekg(pos);
                                break;
                            }else{
                                newInData.setGene(newInData.getGene() + s);
                            }
                        //Get the position of the current(already read) line and keep as to reposition as necessary
                        pos = infile.tellg();
                    }
                    //Add the GeneData to the inputData set once a new FASTA description is discovered and recognized
                    if (newInData.getGene().length() > 0) {
                        inputData.push_back(newInData);
                    }
                }
            }

            //Close and clear input file stream once done
            infile.close();
            infile.clear();
#endif
        //Return the transformed data
        return inputData;
    } else {
        //Throw error in file open
        throw FileInpException();
    }
}

//int File2Obj::getInMSADara(std::string path){
//
//    tesseract::TessBaseAPI tess;
//
//    if (tess.Init("./tessdata", "eng"))
//    {
//        std::cout << "OCRTesseract: Could not initialize tesseract." << std::endl;
//        return 1;
//    }
//
//    // setup
//    tess.SetPageSegMode(tesseract::PageSegMode::PSM_AUTO);
//    tess.SetVariable("save_best_choices", "T");
//
//    // read image
//    auto pixs = pixRead(path.c_str());
//    if (!pixs)
//    {
//        std::cout << "Cannot open input file: " << path.c_str() << std::endl;
//        return 1;
//    }
//
//    // recognize
//    tess.SetImage(pixs);
//    tess.Recognize(0);
//
//    // get result and delete[] returned char* string
//    std::cout << std::unique_ptr<char[]>(tess.GetUTF8Text()).get() << std::endl;
//
//    // cleanup
//    tess.Clear();
//    pixDestroy(&pixs);
//
//    return 0;
//}
