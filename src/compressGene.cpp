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
        case ' ': return 0;
        case 'A': return 4;
        case 'B': return 7;
        case 'C': return 8;
        case 'D': return 9;
        case 'E': return 11;
        case 'F': return 12;
        case 'G': return 14;
        case 'H': return 15;
        case 'I': return 16;
        case 'J': return 17;
        case 'K': return 18;
        case 'L': return 19;
        case 'M': return 22;
        case 'N': return 23;
        case 'O': return 24;
        case 'P': return 25;
        case 'Q': return 27;
        case 'R': return 28;
        case 'S': return 29;
        case 'T': return 30;
        case 'U': return 31;
        case 'V': return 32;
        case 'W': return 34;
        case 'X': return 35;
        case 'Y': return 36;
        case 'Z': return 37;
    }
}

//Return the relevant char to the input int parameter.
//Should be the reverse of the toEncValue conversion.
char toDecValue(int c){
    switch (c){
        case 0: return ' ';
        case 4: return 'A';
        case 7: return 'B';
        case 8: return 'C';
        case 9: return 'D';
        case 11: return 'E';
        case 12: return 'F';
        case 14: return 'G';
        case 15: return 'H';
        case 16: return 'I';
        case 17: return 'J';
        case 18: return 'K';
        case 19: return 'L';
        case 22: return 'M';
        case 23: return 'N';
        case 24: return 'O';
        case 25: return 'P';
        case 27: return 'Q';
        case 28: return 'R';
        case 29: return 'S';
        case 30: return 'T';
        case 31: return 'U';
        case 32: return 'V';
        case 34: return 'W';
        case 35: return 'X';
        case 36: return 'Y';
        case 37: return 'Z';
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
        tempStore+=std::bitset< 6 >( toEncValue(c)).to_string();
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
    tempStore = tempStore.substr(0,tempStore.length()-(tempStore.length()%6));

    //Variable to store result.
    std::string res="";

    //Convert the bit-stream to the original format.
    //Break bit-stream into groups of 6.
    for(int i=0;i<tempStore.length();i+=6){
        //Get a group of 6 bits.
        int j = std::stoi(tempStore.substr(i,6), nullptr, 2);
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
