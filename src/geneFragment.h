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
    std::vector<GeneData> multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex);
    std::vector<GeneData> multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString);
    GeneData indexFragment(GeneData inpData, int strtIndex, int endIndex);
    GeneData stringFragment(GeneData inpData, std::string strtString, std::string endString);
}

#endif //MORFICONVERT_GENEFRAGMENT_H
