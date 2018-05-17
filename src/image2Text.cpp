//
// Created by Damitha on 4/14/2018.
//

#include "image2Text.h"

#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <fstream>
#include <iterator>

#include "GeneData.h"

//Read a bmp file (via input path given) and get the rgb values
std::vector<std::vector<std::vector<char>>> rgbMap(std::string file) {
    //Output data
    std::vector<std::vector<std::vector<char>>> dataOut;

    //Header size of a bmp file
    static constexpr size_t HEADER_SIZE = 54;

    //Open the bmp file
    std::ifstream bmp(file, std::ios::binary);

    //Get the header data
    std::array<char, HEADER_SIZE> header;
    bmp.read(header.data(), header.size());

    //Extract the different header data
    auto fileSize = *reinterpret_cast<uint32_t *>(&header[2]);
    auto dataOffset = *reinterpret_cast<uint32_t *>(&header[10]);
    auto width = *reinterpret_cast<uint32_t *>(&header[18]);
    auto height = *reinterpret_cast<uint32_t *>(&header[22]);
    auto depth = *reinterpret_cast<uint16_t *>(&header[28]);

    //Get the actual bmp image data
    std::vector<char> img(dataOffset - HEADER_SIZE);
    bmp.read(img.data(), img.size());

    //Get the bytes a pixel is represented by
    int pixSet = depth/8;

    //Calculate the total data size
    auto dataSize = ((width * pixSet + pixSet) & (~pixSet)) * height;
    //Resize the img vector via the data size and read the content
    img.resize(dataSize);
    bmp.read(img.data(), img.size());

    for (int i = 0; i < height; i++) {
        //Pixel line
        std::vector<std::vector<char>> dataIOut;
        for (int j = 0; j < width * pixSet; j += pixSet) {
            //A pixel
            std::vector<char> dataJOut;

            //Add the rgb values to the pixel
            dataJOut.push_back(img[j + 2 + i * width* pixSet]);
            dataJOut.push_back(img[j + 1 + i * width* pixSet]);
            dataJOut.push_back(img[j + i * width* pixSet]);
            dataIOut.push_back(dataJOut);
        }
        //As the data in the image is in the reverse order add the data to the start of the ouput vector
        dataOut.insert(dataOut.begin(),dataIOut);
    }

    return dataOut;
}

//Check if given two pixels have the same colour
bool checkIfEqualPix(std::vector<char> pixA, std::vector<char> pixB) {
    //Check all r,g and b values
    if (pixA.at(0) == pixB.at(0)) {
        if (pixA.at(1) == pixB.at(1)) {
            if (pixA.at(2) == pixB.at(2)) {
                //If all are matching return true
                return true;
            }
        }
    }
    //If not matching return false
    return false;
}

//Convert a given "character"'s pixels into a boolean vector where each pixel would either match the background or not
std::vector<bool> getGScale(std::vector<std::vector<char>> cPix){
    //Output
    std::vector<bool> outPixOrien;

    //For the left topmost pixel which should always be the background colour.
    outPixOrien.push_back(false);
    //Iterate through all the pixels in the character
    if(checkIfEqualPix(cPix.at(0),cPix.at(cPix.size()-1))){
        for(int i=1; i<cPix.size(); i++){
            //If the pixel colour matches the background colour pass true, else false
            if(checkIfEqualPix(cPix.at(0),cPix.at(i))){
                outPixOrien.push_back(false);
            } else{
                outPixOrien.push_back(true);
            }
        }
    }
    return outPixOrien;
}

//Get the height and the width of a character in the image
//Takes the number of genes in the image as well as the image itself as parameters
//*Assume that all the heights and the widths of the characters are the same
std::vector<int> getCharMeasure(int geneNum, std::vector<std::vector<std::vector<char>>> inp){
    //Vector that will hold the output ( will only have 2 elements )
    std::vector<int> outCharSize;

    //Calculate height of a character
    int geneHeight = inp.size()/geneNum;

    //Add the height to the output
    outCharSize.push_back(geneHeight);

    //Set the minimum distance at which a colour change occurs(horizontally)
    int minWidth=inp.at(0).size();

    //Iterate through character lines' first line
    for (int i=0; i<inp.size(); i+=geneHeight){
        std::vector<std::vector<char>> imgStrip = inp.at(i);

        //The length at which a change in colour occurs
        int changeLen = 1;
        //Iterate through the line comparing adjacent pixels
        for (int j=1; j<imgStrip.size(); j++){
            //If there is no change then keep on adding to the length. If there is a change then see if it is smaller
            //than the current minWidth and add if so assign the new value else keep the value as it is.
            if(checkIfEqualPix(imgStrip.at(j),imgStrip.at(j-1))){
                changeLen+=1;
            } else {
                minWidth = std::min(minWidth,changeLen);
                changeLen=1;
            }
        }
    }
    //Add the width to the output
    outCharSize.push_back(minWidth);

    return  outCharSize;
}

//Break an image into letter lines / gene and get background colour match for each letter
//Takes the number of genes in the image as well as the image itself as parameters
//*Assume that all the heights and the widths of the characters are the same
std::vector<std::vector<std::vector<bool>>> lineSplit(int geneNum, std::string filePath){
    //Get the image as rbg values
    std::vector<std::vector<std::vector<char>>> imgCharVecs = rgbMap(filePath);

    //Get the measurements of a character
    std::vector<int> getLetterMeasures = getCharMeasure(geneNum, imgCharVecs);

    //Output
    std::vector<std::vector<std::vector<bool>>> outRes;

    //Iterate through each line
    for(int i1=0; i1<imgCharVecs.size(); i1+=getLetterMeasures.at(0)){
        //Resultant letter lines
        std::vector<std::vector<bool>> letterLine;
        //Iterate through the letters
        for(int i2=0; i2<imgCharVecs.at(0).size(); i2+=getLetterMeasures.at(1)){
            //Get the total pixels in one letter
            std::vector<std::vector<char>> letterPic;
            for(int j1=i1;j1<(i1+getLetterMeasures.at(0)); j1++){
                for(int j2=i2;j2<(i2+getLetterMeasures.at(1)); j2++){
                    letterPic.push_back(imgCharVecs.at(j1).at(j2));
                }
            }
            //Get the backgroung colour match and add to the letter line
            letterLine.push_back(getGScale(letterPic));
        }
        //Add to the output
        outRes.push_back(letterLine);
    }

    return outRes;
}

//Convert a vector of boolean values to a string of 1s and 0s.
std::string boolVecToStr(std::vector<bool> eqBoolVec){
    //Output
    std::string equi0n1="";

    for(bool b:eqBoolVec){
        if(b){
            equi0n1+="1";
        } else{
            equi0n1+="0";
        }
    }

    return equi0n1;
}

//Convert an MSA image + gene descriptions to GeneData objects
std::vector<GeneData> Image2Text::getInImgMSAData(std::vector<std::string> genes, std::string filePath){
    //Output
    std::vector<GeneData> outData;

    //Get genes as separate lines (in boolean format)
    std::vector<std::vector<std::vector<bool>>> letterLines = lineSplit(genes.size(),filePath);

    //Iterate through genes
    for (int i=0; i<genes.size();i++){
        std::string geneHeader = genes.at(i);
        std::vector<std::vector<bool>> letterLine = letterLines.at(i);

        //The identified /converted gene sequence
        std::string convGeneSeq="";

        //Iterate through letters
        for(int j=0; j<letterLine.size(); j++){
            std::vector<bool> eqBoolVec = letterLine.at(j);

            //Boolean string of the background colour match
            std::string equi0n1=boolVecToStr(eqBoolVec);

            //Input file stream
            std::ifstream infile;

            //Open file that hold the trained data
            infile.open("Character patterns.morfi");

            //If there was no error in opening the file
            if (infile.good()) {
                //Temporary storage string
                std::string sTemp;

                while (getline(infile, sTemp)) {
                    //Check if the 0 & 1 stings match
                    if (sTemp.substr(sTemp.find_first_of(" ")+1) == equi0n1) {
                        //If they do then add the character represented by them to the gene sequence
                        //*Characters are stored by their ASCII value
                        char addChar = (char)atoi(sTemp.substr(0,sTemp.find_first_of(" ")).c_str());
                        if(addChar!='-'){
                            convGeneSeq+=addChar;
                        }
                        break;
                    }
                }
            }
        }
        //Create a new GeneData object and set its description and gene sequence
        GeneData newInData("", "");
        newInData.setGene(convGeneSeq);
        newInData.setDetails(geneHeader);

        outData.push_back(newInData);
    }

    return outData;
}

//Delete all the trained data
void Image2Text::removeTrainData(){
    remove("Character patterns.morfi");
}

//Train with an MSA image + gene sequences to identify particular letters
void Image2Text::trainImgMSAData(std::vector<std::string> genes, std::string filePath){
    //Get genes as separate lines (in boolean format)
    std::vector<std::vector<std::vector<bool>>> letterLines = lineSplit(genes.size(),filePath);

    //Variable to keep the letters already trained
    std::vector<std::string> addedLetters;

    //Iterate through genes
    for (int i=0; i<genes.size();i++){
        std::string geneSeq = genes.at(i);
        std::vector<std::vector<bool>> letterLine = letterLines.at(i);

        //Iterate through letters
        for(int j=0; j<letterLine.size(); j++){
            //The current letter being trained
            std::string teachNLetter = std::to_string(geneSeq[j]);
            //Boolean values oof the background colour match
            std::vector<bool> eqBoolVec = letterLine.at(j);

            //Boolean string of the background colour match
            std::string equi0n1=boolVecToStr(eqBoolVec);

            //Check if the letter already trained
            bool haveMatch = false;
            for (std::string c:addedLetters){
                if(teachNLetter == c){
                    haveMatch = true;
                    break;
                }
            }

            //Input file stream
            std::ifstream infile;

            //Open file that hold the trained data
            infile.open("Character patterns.morfi");

            //If there was no error in opening the file and untrained
            if (infile.good() && !(haveMatch)) {
                //Temporary storage string
                std::string sTemp;

                //Check if there is a match in the trained data
                while (getline(infile, sTemp)) {
                    if (sTemp.substr(0,sTemp.find_first_of(" ")) == teachNLetter && sTemp.substr(sTemp.find_first_of(" ")+1) == equi0n1) {
                        haveMatch = true;
                        break;
                    }
                }
            }

            //If the data being entered is untrained then add the data to the trained data file
            if(!haveMatch){
                std::ofstream outfile;
                outfile.open("Character patterns.morfi",  std::ios::out | std::ios::app);
                outfile << teachNLetter + " " + equi0n1 + "\n";
                addedLetters.push_back(teachNLetter);
            }

        }
    }
}
