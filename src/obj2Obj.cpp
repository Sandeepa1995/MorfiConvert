//
// Created by Damitha on 3/10/2018.
//

#include "obj2Obj.h"

#include <iostream>
#include <vector>
#include <omp.h>

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

//Incorrect input type exception
class InvalidOutException : public MorfiObj2ObjException {
public:
    virtual const char *what() const throw() {
        return "Invalid data to the output data format type";
    }
};

//Return the data as a set of Vectors
std::vector<std::string> Obj2Obj::giveOutData(std::vector<GeneData> res, std::string outType, int thrdCnt) {
    //Initialize a output vector set
    std::vector<std::string> outVec;

    //**Check if output data type is FASTA**
    if (outType == "FASTA") {
        //Iterate through the GeneData set and write the data to the output vector set
        //OpenMP line to read the data content in GeneData object format using "thrdCnt" number of threads.
#pragma  omp parallel for num_threads(thrdCnt)
        for (int x = 0; x < res.size(); x++) {
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
            //Make adding the new strings to the output vector thread safe.
#pragma omp critical
            {
                outVec.push_back(sTemp);
            }
        }
        return outVec;

        //**Check if the output data type is the list of NCBI identifiers**
    } else if (outType == "NCBI-IDs") {
        //Iterate through the list of gene sequences
        //OpenMP line to read the data content in GeneData object format using "thrdCnt" number of threads.
#pragma  omp parallel for num_threads(thrdCnt)
        for (int x = 0; x < res.size(); x++) {
            GeneData value = res.at(x);
            //Temporary variable to hold the gene description
            std::string gene = value.getDetails();
            //Make adding the new strings to the output vector thread safe.
#pragma omp critical
            {
                outVec.push_back(SupMorfi::getIdentifier(gene));
            }
        }
        return outVec;

        //**Check if output data type is GCG**
    } else if (outType == "GCG") {
        //iterate through the GeneData set and write the data to the output vector set
#pragma  omp parallel for num_threads(thrdCnt)
        for (int x = 0; x < res.size(); x++) {
            //Variables to hold the header and the gene type
            //ex - "!!NA_SEQUENCE 1.0"
            std::string seqHeader;
            //ex - "N"
            std::string geneType;

            GeneData value = res.at(x);

            //Set the header and the gene type
            if (SupMorfi::checkIsNucleic(value)) {
                seqHeader = "!!NA_SEQUENCE 1.0";
                geneType = "N";
            } else {
                seqHeader = "!!AA_SEQUENCE 1.0";
                geneType = "P";
            }

            //Temporary string to store the result that is being built
            std::string sTemp = "";

            //Add header
            sTemp += seqHeader;
            sTemp += "\n\n";

            //Get the identifier and the reset of the description separated from the GeneData object
            std::vector<std::string> splitDes = SupMorfi::separateIdentifier(value.getDetails());

            //Add the description other than the ID
            sTemp += splitDes.at(1);
            sTemp += "\n\n";

            //Add the ID and gene summary
            sTemp += (splitDes.at(0) + "  Length: " + std::to_string(value.getGene().length()) + "  Type: " + geneType +
                      "  Check: " + std::to_string(SupMorfi::gcgChecksum(value.getGene())) + " ..");

            //Variable to hold the i th character of the gene written
            int i = 0;
            for (char c : value.getGene()) {
                //Separate to a new line every 50 characters
                if (i % 50 == 0) {
                    sTemp += "\n\n";
                    i += 1;

                    //Pad the line so that the numbers are aligned left
                    int padding = std::to_string(value.getGene().length()).length() - std::to_string(i).length();
                    for (int k = 0; k <= padding; k++) {
                        sTemp += " ";
                    }

                    //Add the index whithin the gene sequence of the current character
                    sTemp += std::to_string(i);
                    sTemp += " ";
                } else if (i % 10 == 0) {
                    //Separate with a space every 10th character
                    i += 1;
                    sTemp += " ";
                } else {
                    i += 1;
                }
                //Add the character
                sTemp += c;
            }
            sTemp += "\n";
            //Make adding the new strings to the output vector thread safe.
#pragma omp critical
            {
                outVec.push_back(sTemp);
            }
        }
        return outVec;


        //**Check if output data type is EMBL**
    } else if (outType == "EMBL") {
        //EMBL is only valid for nucleic sequences
        if(!SupMorfi::checkAllNucleic(res)){
            throw InvalidOutException();
        }

        //iterate through the GeneData set and write the data to the output vector set
        //OpenMP line to read the data content in GeneData object format using "thrdCnt" number of threads.
#pragma  omp parallel for num_threads(thrdCnt)
        for (int x = 0; x < res.size(); x++) {
            GeneData value = res.at(x);

            //Get the identifier and the reset of the description separated from the GeneData object
            std::vector<std::string> splitDes = SupMorfi::separateIdentifier(value.getDetails());

            //Temporary string to store the result that is being built
            std::string sTemp = ("ID " + splitDes.at(0) + "; SV 1; linear; unassigned DNA; STD; UNC; " +
                                 std::to_string(value.getGene().length()) + " BP.");

            //Get the count of each Nucleic acid
            std::vector<int> countNeuc = SupMorfi::countNeuc(value);
            sTemp += "\n";
            sTemp += "XX";
            sTemp += "\n";

            //Add the description other than the ID
            sTemp += ("DE " + splitDes.at(1));
            sTemp += "\n";
            sTemp += "XX";
            sTemp += "\n";
            sTemp += "XX";
            sTemp += "\n";

            //Add the nucleic acid count
            sTemp += ("SQ Sequence " + std::to_string(value.getGene().length()) + " BP; " +
                      std::to_string(countNeuc.at(0)) + " A; " + std::to_string(countNeuc.at(1)) + " C; " +
                      std::to_string(countNeuc.at(2)) + " G; " + std::to_string(countNeuc.at(3)) + " T; " +
                      std::to_string(countNeuc.at(4)) + " other;" + "\n");

            //Variable to hold the i th character of the gene written
            int i = 0;
            for (char c : value.getGene()) {
                if (i == (value.getGene().length() - 1)) {
                    //Ending of a gene sequence
                    sTemp += c;
                    i += 1;
                    sTemp += " ";
                    sTemp += std::to_string(i);
                    sTemp += "\n";
                } else if ((i % 60 == 0 && i != 0)) {
                    //Separate to a new line every 60 characters
                    sTemp += " ";
                    sTemp += std::to_string(i);
                    sTemp += "\n ";
                    sTemp += c;
                    i += 1;
                } else if (i % 10 == 0) {
                    //Separate with a space every 10th character
                    i += 1;
                    sTemp += " ";
                    sTemp += c;
                } else {
                    i += 1;
                    sTemp += c;
                }
            }
            //Add symbol of gene ending
            sTemp += "//";
            //Make adding the new strings to the output vector thread safe.
#pragma omp critical
            {
                outVec.push_back(sTemp);
            }
        }
        return outVec;

        //**Check if output data type is Genebank**
    } else if (outType == "GenBank") {
        //Iterate through the GeneData set and write the data to the output vector set
        //OpenMP line to read the data content in GeneData object format using "thrdCnt" number of threads.
#pragma  omp parallel for num_threads(thrdCnt)
        for (int x = 0; x < res.size(); x++) {
            GeneData value = res.at(x);

            //Get the identifier and the reset of the description separated from the GeneData object
            std::vector<std::string> splitDes = SupMorfi::separateIdentifier(value.getDetails());
            //Temporary string to store the result that is being built
            std::string sTemp = ("LOCUS " + splitDes.at(0) + " " +
                                 std::to_string(value.getGene().length()) + " bp DNA linear UNC");
            sTemp += "\n";
            //Add the description other than the ID
            sTemp += ("DEFINITION " + splitDes.at(1) + " .");
            sTemp += "\n";

            //Add the ID
            sTemp += ("ACCESSION " + splitDes.at(0).substr(0,splitDes.at(0).find_first_of('.')));
            sTemp += "\n";

            //Add to show that the gene sequence is starting
            sTemp += ("ORIGIN");

            //Variable to hold the i th character of the gene written
            int i = 0;
            for (char c : value.getGene()) {
                if (i % 60 == 0) {
                    //Separate to a new line every 60 characters
                    sTemp += "\n ";
                    i += 1;
                    sTemp += std::to_string(i);
                    sTemp += " ";
                } else if (i % 10 == 0) {
                    //Separate with a space every 10th character
                    i += 1;
                    sTemp += " ";
                } else {
                    i += 1;
                }
                sTemp += c;
            }
            sTemp += "\n";
            //Add symbol of gene ending
            sTemp += "//";
            //Make adding the new strings to the output vector thread safe.
#pragma omp critical
            {
                outVec.push_back(sTemp);
            }
        }
        return outVec;

        //**Check if output data type is PIR or NBRF**
    } else if (outType == "PIR" ||outType == "NBRF" ) {
        //Iterate through the GeneData set and write the data to the output vector set
        //OpenMP line to read the data content in GeneData object format using "thrdCnt" number of threads.
#pragma  omp parallel for num_threads(thrdCnt)
        for (int x = 0; x < res.size(); x++) {
            //Variables to hold the header and the gene type
            // ex - "D1"
            std::string geneType;

            GeneData value = res.at(x);

            //Set the gene type
            if (SupMorfi::checkIsNucleic(value)) {
                geneType = "D1";
            } else {
                geneType = "P1";
            }

            //Get the identifier and the reset of the description separated from the GeneData object
            std::vector<std::string> splitDes = SupMorfi::separateIdentifier(value.getDetails());

            //Temporary string to store the result that is being built
            std::string sTemp = "";
            //Add the ID
            sTemp += ">" + geneType + ";" + splitDes.at(0);
            sTemp += "\n";

            //Add the description other than the ID
            sTemp += splitDes.at(1);
            //Add the bases / length of the string
            sTemp += " , " + std::to_string(value.getGene().length()) + " bases";

            //Variable to hold the i th character of the gene written
            int i = 0;
            for (char c : value.getGene()) {
                if (i % 50 == 0) {
                    //Separate to a new line every 50 characters
                    sTemp += "\n";
                }
                if (i % 10 == 0) {
                    //Separate with a space every 10th character
                    sTemp += " ";
                }
                i += 1;
                sTemp += c;
            }
            sTemp += "*\n";
            //Make adding the new strings to the output vector thread safe.
#pragma omp critical
            {
                outVec.push_back(sTemp);
            }
        }
        return outVec;

        //**Check if output data type is PHYLIP**
    } else if (outType == "PHYLIP") {
        std::string geneType;

        //Iterate through the GeneData set and write the data to the output vector set
        //OpenMP line to read the data content in GeneData object format using "thrdCnt" number of threads.
#pragma  omp parallel for num_threads(thrdCnt)
        for (int x = 0; x < res.size(); x++) {
            GeneData value = res.at(x);
            //Get the NCBI identifier
            std::string iD = SupMorfi::getIdentifier(value.getDetails());
            int padding = iD.length();

            //Temporary string to store the result that is being built
            std::string sTemp = "";
            sTemp += iD;

            //Variable to hold the i th character of the gene written
            int i = 0;
            for (char c : value.getGene()) {
                if (i % 50 == 0 && i!=0) {
                    //Separate to a new line every 50 characters
                    sTemp += "\n";
                    for(int k=0; k<padding;k++){
                        sTemp += " ";
                    }
                } else if (i % 10 == 0 && i!=0) {
                    //Separate with a space every 10th character
                    sTemp += " ";
                }
                i += 1;
                sTemp += c;
            }
            //Make adding the new strings to the output vector thread safe.
#pragma omp critical
            {
                outVec.push_back(sTemp);
            }
        }
        return outVec;
        //Check if the output data type is the list of NCBI identifiers
    } else {
        //Throw error for unknown type
        throw UnknownOutException();
    }
}
