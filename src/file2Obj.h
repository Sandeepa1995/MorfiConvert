#ifndef MORFICONVERT_FILE2OBJ_H
#define MORFICONVERT_FILE2OBJ_H

#include "GeneData.h"
#include <vector>
#include <iostream>

namespace File2Obj {
    std::vector<GeneData> getInNCBIBLASTData(std::string path);
    std::vector<GeneData> getInCustomBLASTData(std::string path);
    std::vector<GeneData> getInCustomMSAData(std::string path);
    std::vector<GeneData> getInFASTAData(std::string path);
//    int getInMSADara(std::string path);

}

#endif