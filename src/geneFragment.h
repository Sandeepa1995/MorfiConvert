//
// Created by Damitha on 3/19/2018.
//

#ifndef MORFICONVERT_GENEFRAGMENT_H
#define MORFICONVERT_GENEFRAGMENT_H

#include <iostream>
#include <vector>
#include <string>

#include "GeneData.h"

namespace GeneFragment{
    //Fragment the given genes' gene sequence by index
    std::vector<GeneData> multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex);

    //Fragment the given genes' gene sequence by strings
    std::vector<GeneData> multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString);

    //Fragment the given gene's gene sequence by index
    GeneData indexFragment(GeneData inpData, int strtIndex, int endIndex);

    //Fragment the given gene's gene sequence by strings
    GeneData stringFragment(GeneData inpData, std::string strtString, std::string endString);
}

#endif //MORFICONVERT_GENEFRAGMENT_H
