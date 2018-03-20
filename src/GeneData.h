//
// Created by Damitha on 3/4/2018.
//

#ifndef MORFICONVERT_GENEDATA_H
#define MORFICONVERT_GENEDATA_H

#include <iostream>

class GeneData {
    std::string gene, details;
public:
    GeneData(std::string, std::string);

    //Setters
    void setGene(std::string);

    void setDetails(std::string);

    //Getters
    std::string getGene();

    std::string getDetails();
};

#endif //MORFICONVERT_GENEDATA_H
