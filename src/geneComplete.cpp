//
// Created by Damitha on 3/11/2018.
//

#include "geneComplete.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <exception>
#include <omp.h>

#include "supplimentMorfi.h"
#include "compressGene.h"

// Exception (super) class for MorfiConvert File2Obj
class MorfiGeneCompleteException : public std::exception {
public:
    virtual const char *what() const throw() {
        return "Error when executing local database functions.";
    }
};

//Input file opening custom exception
class FileInpException : public MorfiGeneCompleteException {
public:
    virtual const char *what() const throw() {
        return "Error in opening input file. Please check if the given path is correct.";
    }
};

//Incorrect input file type custom exception
class FileTypeException : public MorfiGeneCompleteException {
public:
    virtual const char *what() const throw() {
        return "Error in the file type given. Did not contain data in the FASTA format.";
    }
};

//The sorting function when locally saving input gene sequences.
//Take two gene sequence descriptions as parameters
bool sortGene(std::string x, std::string y) {
    //Sort by the NCBI identifier in the string / gene sequence's description.
    return (SupMorfi::getIdentifier(x) < SupMorfi::getIdentifier(y));
}

//Delete the existing local data file
void GeneComplete::deleteExistingCopy(std::string path) {
    //Temporary string to store the records in "points.morfi"
    std::string s;

    //Remove a duplicate file if exists
    remove(path.c_str());

    //To open "points.morfi"
    std::ifstream infile;
    infile.open("points.morfi");

    //Hold valid lines. ( All others than the one deleted )
    std::vector<std::string> validHolder;

    if (infile.good()) {
        while (getline(infile, s)) {
            //Get the part of the record without the pointer location of the end of the file
            std::string shaveLastCheck = s.substr(0, s.find_last_of(" "));
            //Gte the path of the local data file
            std::string dbPathCheck = shaveLastCheck.substr(0, shaveLastCheck.find_last_of(" "));

            //Add the record to the valid set so long as it does not match the path stated
            if (dbPathCheck != path) {
                validHolder.push_back(s);
            }
        }
        infile.clear();
        infile.close();
        //Remove the current "points.morfi" file
        remove("points.morfi");

        //Write the valid records to a new "points.morfi" file
        std::ofstream outPoints;
        outPoints.open("points.morfi", std::ofstream::out | std::ofstream::app);

        for (std::string sPoints:validHolder) {
            outPoints << sPoints + "\n";
        }

        outPoints.close();
    }
}

//Locally save the data in the input FASTA file in format which both compact and can be easily accessed.
void GeneComplete::configureFile(std::string path, std::string outPath) {
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    //Content after the last '\\' character in the input path.
    std::string lastDash = path.substr(path.find_last_of("\\") + 1);
    //Get the 'name' of the database file.
    std::string fileName = lastDash.substr(0, lastDash.find_last_of("."));

    //Position pointer to traverse the file.
    std::streampos pos = 0;
    //To get the position of the end of the file.
    std::streampos maxPnt = 0;

    //If there was no error in opening the file
    if (infile.good()) {
        //Temporary storage string
        std::string s;

        //Count the number of genes stored.
        long i = 0;
        //To count the number of segments the original FASTA database file is broken into.
        int j = 0;
        //To count the maximum length of a gene sequence in a .morfi local database file
        int maxGene = 0;
        //Temporary variable to store the gene sequences and sort them.
        std::vector<std::string> geneVec;

        //Iterate through the database file
        while (getline(infile, s)) {
            //Check for the start of a FASTA description
            if (s[0] == '>') {
                //Get the identifier (NCBI etc.)  from the description
                std::string ncbiID = SupMorfi::getIdentifier(s.substr(1));

                //Variable to hold the gene sequence
                std::string gene;
                //Iterate through ther file and get the gene sequence till the begining of a new sequence is found
                while (getline(infile, s)) {
                    //If a start of another gene sequence (the description) is found reposition the pointer and break through the loop
                    if (s[0] == '>') {
                        infile.clear();
                        infile.seekg(pos);
                        break;
                    } else {
                        gene += s;
                    }
                    //Save the current position within the file for repositioning
                    pos = infile.tellg();
                }

                //Hold the value of the longest gene sequence in the file to file (Needed when retrieving the stored data)
                if (gene.length() > maxGene) {
                    maxGene = gene.length();
                }
                //Push the recognized data onto a vector
                geneVec.push_back(ncbiID + " " + gene);
                //Increase the counter for the gene number
                i++;
            }
            //If 5,000,000 genes have been found write a file to hold the current values.
            //Breaking is done at this point to reduce the time spent sorting and to ensure that the memory used would be in a acceptable range
            if (i == 5000000) {
                //Reset gene counter
                i = 0;
                //Iterate the file number for the database counter
                j += 1;
                //Compress the gene sequences while using OpenMP to introduce parallelism
#pragma omp parallel for
                for (int k = 0; k < geneVec.size(); k++) {
                    //Get the identifier of the gene sequence
                    std::string geneHead = geneVec.at(k).substr(0, geneVec.at(k).find_first_of(" "));
                    //Compress the gene sequence
                    std::string geneSeq = CompressGene::encodeGene(
                            geneVec.at(k).substr(geneVec.at(k).find_first_of(" ") + 1));
                    //Replace the current entry
                    geneVec.at(k) = geneHead + " " + geneSeq;
                }
                //Sort the vector containing the gene sequences. Done to make binary search possible for retrieving.
                std::sort(geneVec.begin(), geneVec.end(), sortGene);

                //Initiate the file name for the file
                std::string fname;
                //Check if a custom output path has been given
                if (outPath != "") {
                    //Build the output file path if custom
                    fname = outPath + "\\" + fileName + std::to_string(j) + ".morfi";
                } else {
                    //Build the output file path if not custom
                    fname = fileName + std::to_string(j) + ".morfi";
                }

                //Delete the record of the if it exists in "points.morfi" and the file itself
                deleteExistingCopy(fname);

                //Open the output file and write the data in the vector
                std::ofstream outfile;
                outfile.open(fname, std::ofstream::out | std::ofstream::trunc);
                for (std::string recrd:geneVec) {
                    outfile << recrd + "\n";
                    //To hold the latest pointer. To be used when retrieving the data
                    maxPnt = outfile.tellp();
                }

                //Output file stream
                std::ofstream outPoints;
                //Open 'points.morfi' file which will keep the metadata of each .morfi local database file
                outPoints.open("points.morfi", std::ofstream::out | std::ofstream::app);

                //Write the metadata of the current file
                outPoints
                        << fname + " " + std::to_string(maxGene) + " " + std::to_string(maxPnt) +
                           "\n";

                outPoints.close();

                //Reset the stored genes
                geneVec = {};
            }
        }
        //Same function as when having 5,000,000 data inputs
        if (i != 0) {
            //Iterate the file number for the database counter
            j += 1;
            //Compress the gene sequences while using OpenMP to introduce parallelism
#pragma omp parallel for
            for (int k = 0; k < geneVec.size(); k++) {
                //Get the identifier of the gene sequence
                std::string geneHead = geneVec.at(k).substr(0, geneVec.at(k).find_first_of(" "));
                //Compress the gene sequence
                std::string geneSeq = CompressGene::encodeGene(
                        geneVec.at(k).substr(geneVec.at(k).find_first_of(" ") + 1));
                //Replace the current entry
                geneVec.at(k) = geneHead + " " + geneSeq;
            }
            //Sort the vector containing the gene sequences. Done to make binary search possible for retrieving.
            std::sort(geneVec.begin(), geneVec.end(), sortGene);

            //Initiate the file name for the file
            std::string fname;
            //Check if a custom output path has been given
            if (outPath != "") {
                //Build the output file path if custom
                fname = outPath + "\\" + fileName + std::to_string(j) + ".morfi";
            } else {
                //Build the output file path if not custom
                fname = fileName + std::to_string(j) + ".morfi";
            }

            //Delete the record of the if it exists in "points.morfi" and the file itself
            deleteExistingCopy(fname);

            //Open the output file and write the data in the vector
            std::ofstream outfile;
            outfile.open(fname, std::ofstream::out | std::ofstream::trunc);
            for (std::string recrd:geneVec) {
                outfile << recrd + "\n";
                //To hold the latest pointer. To be used when retrieving the data
                maxPnt = outfile.tellp();
            }

            //Output file stream
            std::ofstream outPoints;
            //Open 'points.morfi' file which will keep the metadata of each .morfi local database file
            outPoints.open("points.morfi", std::ofstream::out | std::ofstream::app);

            //Write the metadata of the current file
            outPoints
                    << fname + " " + std::to_string(maxGene) + " " + std::to_string(maxPnt) +
                       "\n";

            outPoints.close();
        }

        //Close and clear input file stream once done
        infile.close();
        infile.clear();

        //If proper data has not been found throw an error.
        //Assume this situation would arise when the incorrect file is given as an input
        if (i == 0 && j == 0) {
            //Throw error
            throw FileTypeException();
        }
    } else {
        //Throw error if there is an error in opening the input file
        throw FileInpException();
    }

}

//Method that is to validate the data in the "points.morfi" file.
void GeneComplete::validateLocalData() {
    //Initiate input file stream
    std::ifstream infile;
    //Temporary string to hold the input read from "points.morfi"
    std::string s;

    //Open "points.morfi"
    infile.open("points.morfi");
    //Hold the valid records in "points.morfi"
    std::vector<std::string> validHolder;
    if (infile.good()) {
        //Read "points.morfi"
        while (getline(infile, s)) {
            //Get the part of the record without the pointer location of the end of the file
            std::string shaveLastCheck = s.substr(0, s.find_last_of(" "));
            //Gte the path of the local data file
            std::string dbPathCheck = shaveLastCheck.substr(0, shaveLastCheck.find_last_of(" "));

            //If the file is existing then add the file record to the accepted records
            std::ifstream chkfile;
            chkfile.open(dbPathCheck);
            if (chkfile.good()) {
                validHolder.push_back(s);
            }
        }
        infile.clear();
        infile.close();
        //Remove the current "points.morfi" file
        remove("points.morfi");

        //Write the valid records to a new "points.morfi" file
        std::ofstream outPoints;
        outPoints.open("points.morfi", std::ofstream::out | std::ofstream::app);

        for (std::string sPoints:validHolder) {
            outPoints << sPoints + "\n";
        }

        outPoints.close();
    }
}

//Function to give the complete form of a gene sequence when given the NCBI ID.
std::string GeneComplete::complete(std::string ncbiID, std::string db) {
    //If the database files to be searched are all or not given an empty string search all database files

    //Input file stream
    std::ifstream infile;

    //Open "points.morfi" that contains the metadata for the local data
    infile.open("points.morfi");

    //Temporary string to read the data in "points.morfi"
    std::string sTemp;

    //To check if the gene sequence has been found.
    bool isFound = false;
    //Variable to store the resultant gene
    std::string geneFound = "";

    //Variable to hold the values read from "point.morfi". Done to allow parallelism.
    std::vector<std::string> localRecords;

    //If there is no error in opening the file
    if (infile.good()) {
        //Read all records in "points.morfi"
        while (getline(infile, sTemp)) {
            localRecords.push_back(sTemp);
        }
        //OpenMP line to allow parallelism
#pragma omp parallel for
        for (int rec = 0; rec < localRecords.size(); rec++) {
            std::string s = localRecords.at(rec);

            //If the database to search is not from all or not given as an empty string
            if (db != "all" && db != "") {
                //Check if the accessed record has the same database name as to the one requested
                if (s.substr(s.find_last_of("\\") + 1, db.size()) == db && !isFound) {

                    //Get the pointer position of the end of the file ( End of the line )
                    std::streampos pos = atol(s.substr(s.find_last_of(" ")).c_str());

                    //String without the final pointer value
                    std::string shaveLastPoint = s.substr(0, s.find_last_of(" "));
                    //Get the maximum gene length in the file ( Between the file path and the file pointer )
                    int maxGene = atoi(shaveLastPoint.substr(shaveLastPoint.find_last_of(" ")).c_str() + 1);

                    //Get the file path of the local database file
                    std::string dbPath = shaveLastPoint.substr(0, shaveLastPoint.find_last_of(" "));

                    //To open the local data file
                    std::ifstream datafile;
                    datafile.open(dbPath);

                    //Start and end position variables for binary search
                    long strt = 0;
                    long end = pos;

                    //Temporary string to read the data in the local data files
                    std::string sGene;

                    //Search via binary search till the point where the the difference between the start and the end
                    //variable is less than twice of the maximum length of a gene in the file
                    while (end - strt > maxGene * 2) {
                        //Get mid point
                        long mid = ceil(strt + end) / 2;
                        //Travel to the position in the file given by the mid pointer
                        datafile.clear();
                        datafile.seekg(mid);
                        //Skip the current line ( As in the middle of it )
                        getline(datafile, sGene);
                        //Read the following full line
                        getline(datafile, sGene);
                        //Binary search comparison
                        if (ncbiID < sGene.substr(0, sGene.find_first_of(" "))) {
                            //If the lexicographical value is less in the ID to be found than the ID of the current line
                            //then make the mid point the next end position
                            end = mid;
                        } else if (ncbiID > sGene.substr(0, sGene.find_first_of(" "))) {
                            //If the lexicographical value is more in the ID to be found than the ID of the current line
                            //then make the mid point the next start position
                            strt = mid;
                        } else {
                            //If a match is found then return the value
                            //Prevent collisions from shared access by marking as a critical region
#pragma omp critical
                            {
                                //If the gene found is longer than the existing found gene
                                //Done as files are check in parallel
                                if (geneFound.length() < sGene.substr(sGene.find_first_of(" ") + 1).length()) {
                                    geneFound = sGene.substr(sGene.find_first_of(" ") + 1);
                                }
                                //OpenMP applies a thread pool. Thus by keeping a variable to track if the full gene
                                // has been found, opening later files can be avoided.
                                isFound = true;
                            }
                            break;
                        }
                    }
                    //Reposition within the file
                    datafile.clear();
                    datafile.seekg(strt);

                    //Iterate through the file line by line as the range searched has been narrowed
                    while (getline(datafile, sGene) && datafile.tellg() < end) {
                        //If a match is found then return the value
                        if (ncbiID == sGene.substr(0, sGene.find_first_of(" "))) {
//Prevent collisions from shared access by marking as a critical region
#pragma omp critical
                            {
                                //If the gene found is longer than the existing found gene
                                //Done as files are check in parallel
                                if (geneFound.length() < sGene.substr(sGene.find_first_of(" ") + 1).length()) {
                                    geneFound = sGene.substr(sGene.find_first_of(" ") + 1);
                                }
                                //OpenMP applies a thread pool. Thus by keeping a variable to track if the full gene
                                // has been found, opening later files can be avoided.
                                isFound = true;
                            }
                            break;
                        }
                    }
                }
            } else {
                //Get the pointer position of the end of the file ( End of the line )
                std::streampos pos = atol(s.substr(s.find_last_of(" ")).c_str());

                //String without the final pointer value
                std::string shaveLastPoint = s.substr(0, s.find_last_of(" "));
                //Get the maximum gene length in the file ( Between the file path and the file pointer )
                int maxGene = atoi(shaveLastPoint.substr(shaveLastPoint.find_last_of(" ")).c_str() + 1);

                //Get the file path of the local database file
                std::string dbPath = shaveLastPoint.substr(0, shaveLastPoint.find_last_of(" "));

                //To open the local data file
                std::ifstream datafile;
                datafile.open(dbPath);

                //Start and end position variables for binary search
                long strt = 0;
                long end = pos;

                //Temporary string to read the data in the local data files
                std::string sGene;

                //Search via binary search till the point where the the difference between the start and the end
                //variable is less than twice of the maximum length of a gene in the file
                while (end - strt > maxGene * 2) {
                    //Get mid point
                    long mid = ceil(strt + end) / 2;
                    //Travel to the position in the file given by the mid pointer
                    datafile.clear();
                    datafile.seekg(mid);
                    //Skip the current line ( As in the middle of it )
                    getline(datafile, sGene);
                    //Read the following full line
                    getline(datafile, sGene);
                    //Binary search comparison
                    if (ncbiID < sGene.substr(0, sGene.find_first_of(" "))) {
                        //If the lexicographical value is less in the ID to be found than the ID of the current line
                        //then make the mid point the next end position
                        end = mid;
                    } else if (ncbiID > sGene.substr(0, sGene.find_first_of(" "))) {
                        //If the lexicographical value is more in the ID to be found than the ID of the current line
                        //then make the mid point the next start position
                        strt = mid;
                    } else {
                        //If a match is found then return the value
                        //Prevent collisions from shared access by marking as a critical region
#pragma omp critical
                        {
                            //If the gene found is longer than the existing found gene
                            //Done as files are check in parallel
                            if (geneFound.length() < sGene.substr(sGene.find_first_of(" ") + 1).length()) {
                                geneFound = sGene.substr(sGene.find_first_of(" ") + 1);
                            }
                            //OpenMP applies a thread pool. Thus by keeping a variable to track if the full gene
                            // has been found, opening later files can be avoided.
                            isFound = true;
                        }
                        break;
                    }
                }
                //Reposition within the file
                datafile.clear();
                datafile.seekg(strt);

                //Iterate through the file line by line as the range searched has been narrowed
                while (getline(datafile, sGene) && datafile.tellg() < end) {
                    //If a match is found then return the value
                    if (ncbiID == sGene.substr(0, sGene.find_first_of(" "))) {
//Prevent collisions from shared access by marking as a critical region
#pragma omp critical
                        {
                            //If the gene found is longer than the existing found gene
                            //Done as files are check in parallel
                            if (geneFound.length() < sGene.substr(sGene.find_first_of(" ") + 1).length()) {
                                geneFound = sGene.substr(sGene.find_first_of(" ") + 1);
                            }
                            //OpenMP applies a thread pool. Thus by keeping a variable to track if the full gene
                            // has been found, opening later files can be avoided.
                            isFound = true;
                        }
                        break;
                    }
                }
            }
        }
    }
    return CompressGene::decodeGene(geneFound);
}