//
// Created by Damitha on 3/10/2018.
//

#include "supplimentMorfi.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include "include/json/json/json.h"
#include "GeneData.h"

//Get the file extension when the path is given as a string
std::string SupMorfi::getExtension(std::string path) {
    //Get the component after the last "\" due to the possibility of "." being in the names of directories
    std::string lastDash = path.substr(path.find_last_of("\\") + 1);
    //Get the component after the last "."
    return lastDash.substr(lastDash.find_last_of(".") + 1);
}

//Check if a given string contains only alphabetic characters
bool SupMorfi::isAlpha(std::string s) {
    //Iterate through the string. If not an alphabetic character then return false.
    for (char c:s) {
        int cha = (int) c;
        if (!((cha <= (int) 'z' && cha >= (int) 'a') || (cha <= (int) 'Z' && cha >= (int) 'A'))) {
            return false;
        }
    }
    return true;
}

//Get the NCBI identifier from the gene description
std::string SupMorfi::getIdentifier(std::string s) {
    //Iterate through the string to get a component with a "."
    for (int i = 0; i < s.length(); ++i) {
        if (s[i] == '.') {
            //Variables to hold the starting and ending positions
            int pre = i - 1;
            int post = i + 1;
            //Iterate forward till a white space (Number / Version portion of the identifier)
            while (s[post] != ' ' && s[post] != '|' && s[post] != ';' && post < s.length()) {
                post += 1;
            }
            //Iterate backward till white space (Actual ID)
            while (s[pre] != ' ' && s[pre] != '|' && s[post] != ';' && pre >= 0) {
                pre -= 1;
            }
            //Check if the selected region is matching with the NCBI ID format and then return the value
            if (isdigit(s[post - 1]) && isalnum(s[pre + 1])) {
                return (s.substr(pre + 1, post - pre - 1));
            }
        }
    }
    //If an NCBI identifier was not found then search for a similar identifier
    //Iterate through the string to get a component with a "_"
    for (int i = 0; i < s.length(); ++i) {
        if (s[i] == '_') {
            //Variables to hold the starting and ending positions
            int pre = i - 1;
            int post = i + 1;
            //Iterate forward till a white space (Number / Version portion of the identifier)
            while (s[post] != ' ' && post <= s.length()) {
                post += 1;
            }
            //Iterate backward till white space (Actual ID)
            while (s[pre] != ' ' && pre >= 0) {
                pre -= 1;
            }
            //Check if the selected region is matching with a format that can be called as an identifier
            if (isalnum(s[post - 1]) && isalnum(s[pre + 1])) {
                return (s.substr(pre + 1, post - pre - 1));
            }
        }
    }
    //In case no identifier is found.
    return "";
}

//Extract only the alphabetical characters from the input string
std::string SupMorfi::extractAlpha(std::string req) {
    std::string res;
    for (char c:req) {
        int cha = (int) c;
        //If the character being checked is an alphabetic character then add it to the output string.
        if (((cha <= (int) 'z' && cha >= (int) 'a') || (cha <= (int) 'Z' && cha >= (int) 'A'))) {
            res += c;
        }
    }
    return res;
}

//Function to calculate the checksum used in the GCG representation
int SupMorfi::gcgChecksum(std::string req) {
    int cSum = 0;
    int indx = 0;
    for (char c:req) {
        indx += 1;
        cSum += (indx * (int) c);
        if (indx == 57) {
            indx = 0;
        }
    }
    return cSum % 10000;
}

//Check if a given GeneData object is a nucleic acid
bool SupMorfi::checkIsNucleic(GeneData gene) {
    for (char c: gene.getGene()) {
        if (!((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T') || (c == 'U'))) {
            return false;
        }
    }
    return true;
}

//Multi object function for checking if the given genes are nucleic acids or not
bool SupMorfi::checkAllNucleic(std::vector<GeneData> geneSet) {
    for (GeneData gene: geneSet) {
        if (!(SupMorfi::checkIsNucleic(gene))) {
            return false;
        }
    }
    return true;
}

//Separate the identifier and the remaining gene description
std::vector<std::string> SupMorfi::separateIdentifier(std::string s) {
    //Vector to give the output
    std::vector<std::string> outRes;

    //In case no identifier is found.
    outRes.push_back("");
    outRes.push_back("");

    //Iterate through the description till a '.' is found.
    for (int i = 0; i < s.length(); ++i) {
        if (s[i] == '.') {
            //Extract the ID
            //Variables to hold the starting and ending positions
            int pre = i - 1;
            int post = i + 1;
            //Iterate forward till a white space (Number / Version portion of the identifier)
            while (s[post] != ' ' && post < s.length()) {
                post += 1;
            }
            //Iterate backward till white space (Actual ID)
            while (s[pre] != ' ' && pre >= 0) {
                pre -= 1;
            }

            //Check if the selected region is matching with the NCBI ID format and then add to the output vector set
            if (isdigit(s[post - 1]) && isalnum(s[pre + 1])) {
                //ID
                std::string iD = s.substr(pre + 1, post - pre - 1);
                //Remaining description
                std::string tempDes = s.substr(post - pre-1);


                outRes.at(0) = iD;
                if(iD.length() != s.length()) {
                    if(tempDes[0]==' '){
                        tempDes = tempDes.substr(tempDes.find_first_not_of(" "));
                    }
                    if(tempDes[tempDes.length()-1]==' '){
                        tempDes = tempDes.substr(0,tempDes.find_last_of(" "));
                    }
                    outRes.at(1) = tempDes;
                }
                return outRes;
            }
        }
    }
    return outRes;
}

//Count the number of instances a nucleic is present in a gene sequence
std::vector<int> SupMorfi::countNeuc(GeneData gene) {
    //Variables to hold the counted number
    int cA = 0;
    int cC = 0;
    int cG = 0;
    int cT = 0;
    int cO = 0;

    //Iterate through the gene
    for (char c:gene.getGene()) {
        if (c == 'A') {
            cA += 1;
        } else if (c == 'C') {
            cC += 1;
        } else if (c == 'G') {
            cG += 1;
        } else if (c == 'T') {
            cT += 1;
        } else {
            cO += 1;
        }
    }

    //Add the result to the output and return
    std::vector<int> res;
    res.push_back(cA);
    res.push_back(cC);
    res.push_back(cG);
    res.push_back(cT);
    res.push_back(cO);
    return res;
}

//Get the length of the longest gene sequence in a GeneData object set
int SupMorfi::getLongestGeneLength(std::vector<GeneData> geneSet) {
    int mxLen = 0;
    for (GeneData gene:geneSet) {
        if (gene.getGene().length() > mxLen) {
            mxLen = gene.getGene().length();
        }
    }
    return mxLen;
}