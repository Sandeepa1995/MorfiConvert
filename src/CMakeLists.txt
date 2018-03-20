cmake_minimum_required(VERSION 3.9)
project(morficonvert)

# The version number.
set (morficonvert_VERSION_MAJOR 1)
set (morficonvert_VERSION_MINOR 0)

# Set the cpp standard
set(CMAKE_CXX_STANDARD 17)

if (MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

#Pre-made
# Check if OpenMP is supported by the compiler
# If so then assign flags accordingly
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

# Add the Jsoncpp library to read json files easily (Add directories)
include_directories("${CMAKE_SOURCE_DIR}/include/json"
        "${CMAKE_SOURCE_DIR}/include/json/json")
add_subdirectory(include/json)
#
## Add the UnQLite library to store the positions of gene data in particular files (Add directories)
#include_directories("${CMAKE_SOURCE_DIR}/include/unqlite")
#add_subdirectory(include/unqlite)
#
##set(CPPAN_BUILD_SHARED_LIBS ON)
#add_subdirectory(.cppan)

# Set the files in the library
set(LIBRARY_FILES
        morficonvert.cpp
        morficonvert.h
        file2Obj.cpp
        file2Obj.h
        GeneData.cpp
        GeneData.h
        obj2Obj.cpp
        obj2Obj.h
        obj2File.cpp
        obj2File.h
        supplimentMorfi.cpp
        supplimentMorfi.h
        geneComplete.cpp
        geneComplete.h
        geneFragment.cpp
        geneFragment.h)

#Add the files in the library
add_library(morficonvert ${LIBRARY_FILES})

# Pre-made linking
# Link the OpenMP library
target_link_libraries(morficonvert -fopenmp)
target_link_libraries(morficonvert Json)
#target_link_libraries(morficonvert unqlite)
#target_link_libraries(morficonvert
#        pvt.cppan.demo.google.tesseract.libtesseract
#        pvt.cppan.demo.danbloomberg.leptonica
#        )