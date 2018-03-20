//
// Created by Damitha on 3/19/2018.
//

#include "geneFragment.h"

#include <iostream>
#include <vector>
#include <string>

#include "GeneData.h"

std::vector<GeneData> GeneFragment::multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex) {
    std::vector<GeneData> outData;

    for (GeneData gene:inpData) {
        GeneData newTempgene("", "");
        newTempgene.setDetails(gene.getDetails());
        newTempgene.setGene(gene.getGene().substr(strtIndex, endIndex-strtIndex));
        outData.push_back(newTempgene);
    }

    return outData;
}

GeneData GeneFragment::indexFragment(GeneData inpData, int strtIndex, int endIndex) {
    GeneData newTempgene("", "");
    newTempgene.setDetails(inpData.getDetails());
    newTempgene.setGene(inpData.getGene().substr(strtIndex, endIndex-strtIndex));
    return newTempgene;
}


std::vector<GeneData>
GeneFragment::multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString) {
    std::vector<GeneData> outData;

    for (GeneData gene:inpData) {
        GeneData newTempgene("", "");
        newTempgene.setDetails(gene.getDetails());
        newTempgene.setGene(gene.getGene().substr(gene.getGene().find_first_of(strtString),
                                                  gene.getGene().find_last_of(endString)-gene.getGene().find_first_of(strtString)));
        outData.push_back(newTempgene);
    }
    return outData;
}

GeneData
GeneFragment::stringFragment(GeneData inpData, std::string strtString, std::string endString){
    GeneData newTempgene("", "");
    newTempgene.setDetails(inpData.getDetails());
    newTempgene.setGene(inpData.getGene().substr(inpData.getGene().find_first_of(strtString),
                                              inpData.getGene().find_last_of(endString)-inpData.getGene().find_first_of(strtString)));
    return newTempgene;
}