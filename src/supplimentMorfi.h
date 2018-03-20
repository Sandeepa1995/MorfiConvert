//
// Created by Damitha on 3/10/2018.
//

#ifndef MORFICONVERT_SUPPLIMENTMORFI_H
#define MORFICONVERT_SUPPLIMENTMORFI_H

#include <iostream>
#include <string>
#include <vector>

namespace SupMorfi{
    std::string getExtension(std::string path);
    bool isAlpha(std::string s);
    std::string identify(std::string path);
    std::string getIdentifier(std::string s);
}

#endif //MORFICONVERT_SUPPLIMENTMORFI_H
