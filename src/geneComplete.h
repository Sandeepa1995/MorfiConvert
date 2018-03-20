//
// Created by Damitha on 3/11/2018.
//

#ifndef MORFICONVERT_GENECOMPLETE_H
#define MORFICONVERT_GENECOMPLETE_H

#include <iostream>
#include <string>

namespace GeneComplete{
    void configureFile(std::string path);
    std::string complete(std::string ncbiID,std::string path);
}

#endif //MORFICONVERT_GENECOMPLETE_H
