//
// Created by Damitha on 3/2/2018.
//
#include <iostream>
#include "src/morficonvert.h"
#include "src/GeneData.h"
#include <fstream>
#include <sstream>

int main(){
    std::string path = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Valid Total (NCBI)\\Test1.xml";
    std::string cusPath = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Valid (Custom)\\Test1.txt";

    std::string datapath = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Data\\env_nr\\env_nr";
    std::string imgpath = "F:\\Downloads\\tesseract_example-master\\tesseract_example-master\\with_cmake\\img\\phototest.tif";

    std::string wrdpath = "F:\\Achademic\\CS\\Semester 05\\ABD.docx";
    std::string msatpath = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\MSA Text\\Test1.txt";
//    std::vector<GeneData> test= morfiConvert::getInData(path, "NCBI-BLAST");
//    std::vector<GeneData> test= morfiConvert::getInData(cusPath, "BLAST");
    std::vector<GeneData> test= morfiConvert::getInData(msatpath, "MSA");
    for (GeneData i: test) {
        std::cout << ">" << i.getDetails() << std::endl;
        std::cout << i.getGene() << std::endl;
    }
//    std::cout << test.size();
//    morfiConvert::configLocal(datapath);
//    std::cout << morfiConvert::completeGene("ARO40708.1","env_nr") <<std::endl;
//    std::ifstream infile;
//
//    //Open file
//    infile.open(wrdpath);
//    std::string s;
//
//    //If there was no error in opening the file
//    while (infile.good()) {
//        getline(infile,s);
//        std::cout<<s<<std::endl;
//    }

    return 0;

}