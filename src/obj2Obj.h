//
// Created by Damitha on 3/10/2018.
//

#ifndef MORFICONVERT_OBJ2OBJ_H
#define MORFICONVERT_OBJ2OBJ_H

#include <iostream>
#include <vector>

#include "GeneData.h"

namespace Obj2Obj {
    //Return the data as a set of Vectors
    std::vector <std::string> giveOutData(std::vector <GeneData> res, std::string outType, int thrdCnt=10);
}

#endif //MORFICONVERT_OBJ2OBJ_H
