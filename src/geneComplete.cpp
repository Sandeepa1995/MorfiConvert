//
// Created by Damitha on 3/11/2018.
//

#include "geneComplete.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <math.h>

#include "include/unqlite/unqlite.h"
#include "include/json/json/json.h"

#include "supplimentMorfi.h"


bool sortGene(std::string x, std::string y) {
    return (SupMorfi::getIdentifier(x) < SupMorfi::getIdentifier(y));
}

void GeneComplete::configureFile(std::string path) {
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(path);

    std::string lastDash = path.substr(path.find_last_of("\\") + 1);
    std::string fileName = lastDash.substr(0, lastDash.find_last_of("."));

    //Temporary variable to store the position
    std::streampos pos = 0;
    std::streampos maxPnt = 0;

    std::remove(fileName.c_str());

    //If there was no error in opening the file
    if (infile.good()) {
        //Temporary storage string
        std::string s;

        long i = 0;
        int j = 0;
        int maxGene = 0;
        std::vector<std::string> geneVec;

//        std::ifstream copyCurrent;
//        copyCurrent.open("points.morfi");

        std::ofstream outPoints;
        outPoints.open("points.morfi", std::ofstream::out | std::ofstream::app);

        while (getline(infile, s)) {
            if (s[0] == '>') {
                std::string ncbiID = SupMorfi::getIdentifier(s.substr(1));

                std::string gene;
                while (getline(infile, s)) {
                    if (s[0] == '>') {
                        infile.clear();
                        infile.seekg(pos);
                        break;
                    } else {
                        gene += s;
                    }
                    pos = infile.tellg();
                }

                if (gene.length() > maxGene) {
                    maxGene = gene.length();
                }
                geneVec.push_back(ncbiID + " " + gene);
                i++;
            }
            if (i == 5000000) {
                i = 0;
                j += 1;
                std::sort(geneVec.begin(), geneVec.end(), sortGene);

                std::ofstream outfile;
                outfile.open(fileName + std::to_string(j) + ".morfi", std::ofstream::out | std::ofstream::trunc);
                for (std::string recrd:geneVec) {
                    outfile << recrd + "\n";
                    maxPnt = outfile.tellp();
                }
                outPoints
                        << fileName + std::to_string(j) + " " + std::to_string(maxGene) + " " + std::to_string(maxPnt) +
                           "\n";
                geneVec = {};
            }

        }
        if (i != 0) {
            j += 1;
            std::sort(geneVec.begin(), geneVec.end(), sortGene);
            std::ofstream outfile;
            outfile.open(fileName + std::to_string(j) + ".morfi", std::ofstream::out | std::ofstream::trunc);
            for (std::string recrd:geneVec) {
                outfile << recrd + "\n";
                maxPnt = outfile.tellp();
            }
            outPoints << fileName + std::to_string(j) + " " + std::to_string(maxGene) + " " + std::to_string(maxPnt) +
                         "\n";
        }

        //Close and clear input file stream once done
        infile.close();
        infile.clear();
    }

}

//AGS56207.1 1687030547
std::string GeneComplete::complete(std::string ncbiID, std::string db) {
    if (db != "all") {
        //Input file stream
        std::ifstream infile;

        //Open file
        infile.open("points.morfi");

        std::string s;
        std::string sGene;

        int j = 1;

        bool isFound = false;
        std::string geneFound = "";

        if (infile.good()) {
            while (getline(infile, s) && !isFound) {
                if (s.substr(0, db.size()) == db) {
                    std::streampos pos = atol(s.substr(s.find_last_of(" ")).c_str());
                    int maxGene = atoi(s.substr(s.find_first_of(" ") + 1, s.find_last_of(" ")).c_str());
                    std::ifstream datafile;
                    datafile.open(db + std::to_string(j) + ".morfi");

                    long strt = 0;
                    long end = pos;

                    while (end - strt > maxGene * 2) {
                        long mid = ceil(strt + end) / 2;
                        datafile.clear();
                        datafile.seekg(mid);
                        getline(datafile, sGene);
                        getline(datafile, sGene);
                        if (ncbiID < sGene.substr(0, sGene.find_first_of(" "))) {
                            end = mid;
                        } else if (ncbiID > sGene.substr(0, sGene.find_first_of(" "))) {
                            strt = mid;
                        } else {
                            geneFound = sGene.substr(sGene.find_first_of(" ") + 2);
                            isFound = true;
                            break;
                        }
                    }
                    datafile.clear();
                    datafile.seekg(strt);
                    while (getline(datafile, sGene) && datafile.tellg() < end) {
                        if (ncbiID == sGene.substr(0, sGene.find_first_of(" "))) {
                            geneFound = sGene.substr(sGene.find_first_of(" ") + 2);
                            isFound = true;
                            break;
                        }
                    }
                    j += 1;
                }
            }
        }
        return geneFound;
    } else {
        //Input file stream
        std::ifstream infile;

        //Open file
        infile.open("points.morfi");

        std::string s;
        std::string sGene;

        int j = 1;

        bool isFound = false;
        std::string geneFound = "";

        if (infile.good()) {
            while (getline(infile, s) && !isFound) {
                std::streampos pos = atol(s.substr(s.find_last_of(" ")).c_str());
                int maxGene = atoi(s.substr(s.find_first_of(" ") + 1, s.find_last_of(" ")).c_str());
                std::ifstream datafile;
                datafile.open(s.substr(0, s.find_first_of(" ")) + ".morfi");

                long strt = 0;
                long end = pos;

                while (end - strt > maxGene * 2) {
                    long mid = ceil(strt + end) / 2;
                    datafile.clear();
                    datafile.seekg(mid);
                    getline(datafile, sGene);
                    getline(datafile, sGene);
                    if (ncbiID < sGene.substr(0, sGene.find_first_of(" "))) {
                        end = mid;
                    } else if (ncbiID > sGene.substr(0, sGene.find_first_of(" "))) {
                        strt = mid;
                    } else {
                        geneFound = sGene.substr(sGene.find_first_of(" ") + 2);
                        isFound = true;
                        break;
                    }
                }
                datafile.clear();
                datafile.seekg(strt);
                while (getline(datafile, sGene) && datafile.tellg() < end) {
                    if (ncbiID == sGene.substr(0, sGene.find_first_of(" "))) {
                        geneFound = sGene.substr(sGene.find_first_of(" ") + 2);
                        isFound = true;
                        break;
                    }
                }
            }
        }
        return geneFound;
    }
}