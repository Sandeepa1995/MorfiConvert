//
// Created by Damitha on 4/14/2018.
//

#ifndef MORFICONVERT_IMAGE2TEXT_H
#define MORFICONVERT_IMAGE2TEXT_H

#include <iostream>
#include <string>
#include <vector>

#include "GeneData.h"

namespace Image2Text {
    //Convert an MSA image + gene descriptions to GeneData objects
    std::vector<GeneData> getInImgMSAData(std::vector<std::string> genes, std::string filePath);

    //Train with an MSA image + gene sequences to identify particular letters
    void trainImgMSAData(std::vector<std::string> genes, std::string filePath);

    //Delete all the trained data
    void removeTrainData();
}


#endif //MORFICONVERT_IMAGE2TEXT_H
