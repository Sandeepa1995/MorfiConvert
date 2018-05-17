//
// Created by Damitha on 4/14/2018.
//

#ifndef MORFICONVERT_COMPRESSGENE_H
#define MORFICONVERT_COMPRESSGENE_H

#include <iostream>
#include <string>
#include <vector>

namespace CompressGene {
    //Encode / Compress the input gene
    std::string encodeGene(std::string gene);

    //Decode / Decompress the input gene
    std::string decodeGene(std::string gene);
}


#endif //MORFICONVERT_COMPRESSGENE_H
