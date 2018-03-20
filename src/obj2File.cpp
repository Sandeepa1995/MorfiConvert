//
// Created by Damitha on 3/10/2018.
//

#include "obj2File.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include "GeneData.h"
#include "supplimentMorfi.h"


//Custom exceptions

// Exception (super) class for MorfiConvert File2Obj
class MorfiObj2FileException : public std::exception {
public:
    virtual const char *what() const throw() {
        return "Error when converting the set of GeneData objects to the output file.";
    }
};

//Output file opening custom exception
class FileOutException : public MorfiObj2FileException {
public:
    virtual const char *what() const throw() {
        return "Error in opening output file. Please check if the given path is correct.";
    }
};

//Incorrect output type exception
class UnknownOutException : public MorfiObj2FileException {
public:
    virtual const char *what() const throw() {
        return "Unsupported output data format type";
    }
};

//Write the GeneData vector set into a file with the specified format
void Obj2File::writeOutData(std::vector<GeneData> res, std::string path, std::string outType = "FASTA") {
    //Check if output data type is FASTA
    if (outType == "FASTA") {

        //Make the file a .txt file if not given
        if (SupMorfi::getExtension(path) != "txt") {
            path += ".txt";
        }

        //Initialize a output file stream and open the file
        std::ofstream outfile;
        outfile.open(path, std::ofstream::out | std::ofstream::trunc);

        //Check if the path given exists and throw an exception if not
        if (outfile.good()) {
            //iterate through the GeneData set and write the data to the opened file
            for (auto &value:res) {
                outfile << ">";
                outfile << value.getDetails();
                //FASTA should be split into lines of a maximum of 80 characters
                int i = 0;
                for (char &c : value.getGene()) {
                    if (i % 80 == 0) {
                        outfile << "\n";
                    }
                    i += 1;
                    outfile << c;
                }
                outfile << "\n";
            }
        } else {
            throw FileOutException();
        }
    } else {
        //Throw error for unkown type
        throw UnknownOutException();
    }
}