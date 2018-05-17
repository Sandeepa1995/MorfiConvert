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
void Obj2File::writeOutData(std::vector<std::string> res, std::string path, std::string outType = "FASTA") {
    //Check the output data type and add an extension to the file if it does not match the standard.
    if (outType == "FASTA") {
        //Make the file a .fasta file if not given
        if (SupMorfi::getExtension(path) != "fasta") {
            path += ".fasta";
        }
    } else if (outType == "GCG") {
        //Make the file a .gcg file if not given
        if (SupMorfi::getExtension(path) != "gcg") {
            path += ".gcg";
        }
    } else if (outType == "GenBank") {
        //Make the file a .genbank file if not given
        if (SupMorfi::getExtension(path) != "genbank") {
            path += ".genbank";
        }
    } else if (outType == "PIR") {
        //Make the file a .pir file if not given
        if (SupMorfi::getExtension(path) != "pir") {
            path += ".pir";
        }
    } else if (outType == "NBRF") {
        //Make the file a .nbrf file if not given
        if (SupMorfi::getExtension(path) != "nbrf") {
            path += ".nbrf";
        }
    } else if (outType == "PHYLIP") {
        //Make the file a .phylip file if not given
        if (SupMorfi::getExtension(path) != "phylip") {
            path += ".phylip";
        }
    } else {
        //Throw error for unkown type
        throw UnknownOutException();
    }

    //Initialize a output file stream and open the file
    std::ofstream outfile;
    outfile.open(path, std::ofstream::out | std::ofstream::trunc);

    //Check if the path given exists and throw an exception if not
    if (outfile.good()) {
        //iterate through the GeneData set and write the data to the opened file
        for (auto &value:res) {
            outfile << value;
            outfile << "\n";
        }
    } else {
        throw FileOutException();
    }
}