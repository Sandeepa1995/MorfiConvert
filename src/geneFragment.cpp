//
// Created by Damitha on 3/19/2018.
//

#include "geneFragment.h"

#include <iostream>
#include <vector>
#include <string>

#include "GeneData.h"

//Fragment the given genes' gene sequence by index
//Takes the input genes as a GeneData vector set and the starting index as well as the ending index of the gene sequence
std::vector<GeneData> GeneFragment::multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex) {
    //Output genes
    std::vector<GeneData> outData;

    //To make sure the start index given by the input is not negative
    strtIndex = std::max(0, strtIndex);

    //Iterate through the genes
    for (GeneData gene:inpData) {
        //New gene object
        GeneData newTempgene("", "");

        //Check and set the ending index. If it is greater than the length of the gene sequence then set it as the
        //length of the gene sequence
        int endIndexTemp;
        if (endIndex > gene.getGene().length()) {
            endIndexTemp = gene.getGene().length();
        } else {
            endIndexTemp = endIndex;
        }

        //Copy the description of the original GeneData object to the new one
        newTempgene.setDetails(gene.getDetails());
        //Get the needed gene component and set the gene
        newTempgene.setGene(gene.getGene().substr(strtIndex, endIndexTemp - strtIndex));
        //Push the new gene to the output vector
        outData.push_back(newTempgene);
    }

    return outData;
}

//Fragment the given gene's gene sequence by index
//Takes the input gene as a GeneData object and the starting index as well as the ending index of the gene sequence
GeneData GeneFragment::indexFragment(GeneData inpData, int strtIndex, int endIndex) {
    //To make sure the start index given by the input is not negative
    strtIndex = std::max(0, strtIndex);

    //Check and set the ending index. If it is greater than the length of the gene sequence then set it as the
    //length of the gene sequence
    if (endIndex > inpData.getGene().length()) {
        endIndex = inpData.getGene().length();
    }

    //New gene object
    GeneData newTempgene("", "");

    //Copy the description of the original GeneData object to the new one
    newTempgene.setDetails(inpData.getDetails());
    //Get the needed gene component and set the gene
    newTempgene.setGene(inpData.getGene().substr(strtIndex, endIndex - strtIndex));

    return newTempgene;
}

//Fragment the given genes' gene sequence by strings
//Takes the input genes as a GeneData vector set and the string matching the start of the gene sequences (first match)
// as well as the string matching the end of the gene sequences (last match)
std::vector<GeneData>
GeneFragment::multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString) {
    //Output genes
    std::vector<GeneData> outData;

    //Iterate through the genes
    for (GeneData gene:inpData) {
        //New gene object
        GeneData newTempgene("", "");

        //Copy the description of the original GeneData object to the new one
        newTempgene.setDetails(gene.getDetails());
        //Get the needed gene component and set the gene
        newTempgene.setGene(gene.getGene().substr(gene.getGene().find_first_of(strtString),
                                                  gene.getGene().rfind(endString) -
                                                  gene.getGene().find_first_of(strtString) + endString.length() - 1));
        //Push the new gene to the output vector
        outData.push_back(newTempgene);
    }
    return outData;
}

//Fragment the given gene's gene sequence by strings
//Takes the input gene as a GeneData object and the string matching the start of the gene sequence (first match)
// as well as the string matching the end of the gene sequence (last match)
GeneData
GeneFragment::stringFragment(GeneData inpData, std::string strtString, std::string endString) {
    //New gene object
    GeneData newTempgene("", "");

    //Copy the description of the original GeneData object to the new one
    newTempgene.setDetails(inpData.getDetails());
    //Get the needed gene component and set the gene
    newTempgene.setGene(inpData.getGene().substr(inpData.getGene().find_first_of(strtString),
                                                 inpData.getGene().rfind(endString) -
                                                 inpData.getGene().find_first_of(strtString) + endString.length()));
    return newTempgene;
}