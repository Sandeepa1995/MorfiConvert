//
// Created by Damitha on 3/10/2018.
//

#ifndef MORFICONVERT_SUPPLIMENTMORFI_H
#define MORFICONVERT_SUPPLIMENTMORFI_H

#include <iostream>
#include <string>
#include <vector>
#include "GeneData.h"

namespace SupMorfi {
    //Get the file extension when the path is given as a string
    std::string getExtension(std::string path);

    //Check if a given string contains only alphabetic characters
    bool isAlpha(std::string s);

    //Get the NCBI identifier from the gene description
    std::string getIdentifier(std::string s);

    //Extract only the alphabetical characters from the input string
    std::string extractAlpha(std::string req);

    //Function to calculate the checksum used in the GCG representation
    int gcgChecksum(std::string req);

    //Check if a given GeneData object is a nucleic acid
    bool checkIsNucleic(GeneData gene);

    //Multi object function for checking if the given genes are nucleic acids or not
    bool checkAllNucleic(std::vector<GeneData> geneSet);

    //Separate the identifier and the remaining gene description
    std::vector<std::string> separateIdentifier(std::string s);

    //Count the number of instances a nucleic is present in a gene sequence
    std::vector<int> countNeuc(GeneData gene);

    //Get the length of the longest gene sequence in a GeneData object set
    int getLongestGeneLength(std::vector<GeneData> geneSet);
}

#endif //MORFICONVERT_SUPPLIMENTMORFI_H
