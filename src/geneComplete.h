//
// Created by Damitha on 3/11/2018.
//

#ifndef MORFICONVERT_GENECOMPLETE_H
#define MORFICONVERT_GENECOMPLETE_H

#include <iostream>
#include <string>

namespace GeneComplete{
    //Locally save the data in the input FASTA file in format which both compact and can be easily accessed.
    void configureFile(std::string path, std::string outPath);

    //Function to give the complete form of a gene sequence when given the NCBI ID.
    std::string complete(std::string ncbiID,std::string path);

    //Method that is to validate the data in the "points.morfi" file.
    void validateLocalData();

    //Delete the existing local data file
    void deleteExistingCopy(std::string path);
}

#endif //MORFICONVERT_GENECOMPLETE_H
