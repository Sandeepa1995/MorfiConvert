//
// Created by Damitha on 4/14/2018.
//

#include "compressGene.h"

#include <iostream>
#include <string>
#include <vector>
#include <bitset>

//Return the relevant int number to the char given as the input parameter.
int toEncValue(char c){
    switch (c){
        case 'A': return 0;
        case 'C': return 4;
        case 'D': return 7;
        case 'E': return 8;
        case 'F': return 9;
        case 'G': return 11;
        case 'H': return 12;
        case 'I': return 14;
        case 'K': return 15;
        case 'L': return 16;
        case 'M': return 17;
        case 'N': return 18;
        case 'O': return 19;
        case 'P': return 22;
        case 'Q': return 23;
        case 'R': return 24;
        case 'S': return 25;
        case 'T': return 27;
        case 'U': return 28;
        case 'V': return 29;
        case 'W': return 30;
        case 'Y': return 31;
    }
}

//Return the relevant char to the input int parameter.
//Should be the reverse of the toEncValue conversion.
char toDecValue(int c){
    switch (c){
        case 0: return 'A';
        case 4: return 'C';
        case 7: return 'D';
        case 8: return 'E';
        case 9: return 'F';
        case 11: return 'G';
        case 12: return 'H';
        case 14: return 'I';
        case 15: return 'K';
        case 16: return 'L';
        case 17: return 'M';
        case 18: return 'N';
        case 19: return 'O';
        case 22: return 'P';
        case 23: return 'Q';
        case 24: return 'R';
        case 25: return 'S';
        case 27: return 'T';
        case 28: return 'U';
        case 29: return 'V';
        case 30: return 'W';
        case 31: return 'Y';
    }
}

//Encode the given gene to a compressed format of approximately 62.5% of the original.
std::string CompressGene::encodeGene(std::string gene){
    //Take input gene;
    //Add a tailing letter to later identify the actual tail of the original gene when decoding.
    std::string finGene = gene + "Y";

    //Temporary store to keep the bit-stream value as a string.
    std::string tempStore="";

    //Encode characters individually.
    for (char c:finGene){
        tempStore+=std::bitset< 5 >( toEncValue(c)).to_string();
    }

    //Pad the tempStore variable so that bit-stream to byte-stream conversion occurs as expected
    while(tempStore.length()%8!=0){
        tempStore+="0";
    }

    //Variable to hold the output result.
    std::string res="";

    //Gather the bits in groups of 8 to convert to a byte and get the relevant ASCII character.
    for(int i=0;i<tempStore.length();i+=8){
        //Group bits in groups of 8 and get the relevant int value.
        int j = std::stoi(tempStore.substr(i,8), nullptr, 2);
        //Convert the int value to the ASCII char value.
        res+=(char)j;
    }

    //Return result.
    return res;
}


//Decompress / Decode the given compressed / encoded input gene sequence  to its original format.
std::string CompressGene::decodeGene(std::string gene){
    //Temporary variable to store the bit-stream.
    std::string tempStore="";

    //Convert characters written in a bit-stream.
    for (char c:gene){
        tempStore+=std::bitset< 8 >( (int)(c)).to_string();
    }

    //Remove the padded bits.
    tempStore = tempStore.substr(0,tempStore.length()-(tempStore.length()%5));

    //Variable to store result.
    std::string res="";

    //Convert the bit-stream to the original format.
    //Break bit-stream into groups of 5.
    for(int i=0;i<tempStore.length();i+=5){
        //Get a group of 5 bits.
        int j = std::stoi(tempStore.substr(i,5), nullptr, 2);
        //Decode / Decompress the int taken.
        char addC = toDecValue(j);
        //If valid add to the output gene.
        if(addC!=' '){
            res+=addC;
        }
    }

    //Remove the tailing 'Y' character added during encoding and return the result.
    return res.substr(0,res.find_last_of('Y'));
}
