//
// Created by Damitha on 3/10/2018.
//

#ifndef MORFICONVERT_OBJ2FILE_H
#define MORFICONVERT_OBJ2FILE_H

#include "GeneData.h"
#include <vector>
#include <iostream>

namespace Obj2File{
    void writeOutData(std::vector<GeneData> res, std::string path, std::string outType);
}

#endif //MORFICONVERT_OBJ2FILE_H
