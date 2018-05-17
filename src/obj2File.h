//
// Created by Damitha on 3/10/2018.
//

#ifndef MORFICONVERT_OBJ2FILE_H
#define MORFICONVERT_OBJ2FILE_H

#include "GeneData.h"
#include <vector>
#include <iostream>

namespace Obj2File{
    //Write the GeneData vector set into a file with the specified format
    void writeOutData(std::vector<std::string> res, std::string path, std::string outType);
}

#endif //MORFICONVERT_OBJ2FILE_H
