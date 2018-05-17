//
// Created by Damitha on 3/4/2018.
//

#include "GeneData.h"

//Gene Data Constructor
GeneData::GeneData(std::string dtls,std::string gen) {
    gene = gen;
    details = dtls;
}

//Setter for the Gene Sequence
void GeneData::setGene(std::string gen) {
    gene = gen;
}

//Setter for the Gene-data details
void GeneData::setDetails(std::string dtls) {
    details = dtls;
}

//Getter for the Gene Sequence
std::string GeneData::getGene() {
    return gene;
}

//Getter for the Gene-data details
std::string GeneData::getDetails() {
    return details;
}