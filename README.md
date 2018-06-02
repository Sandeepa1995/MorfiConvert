# MorfiConvert : Bioinformatics Data Format Converter

MorfiConvert is a bioinformatics data format converter created to specialize in converting the outputs of one bioinformatics process to a standard format that can be used as an input for another bioinformatics process.

# Using MorfiConvert 
MorfiConvert is uses CMake to build its file to give itself platform independency. It also uses OpenMP to support parallelism. Thus, to use MorfiConvert you will need to,
1.	Have a C++ compiler (that has OpenMP)  and install CMake
2.	Build the project via entering the following in the command-line within the MorfiConvert folder,
```sh
mkdir build && cd build
cmake ..
```
or 

Using the MorfiConvert in a project where CMake is used via downloading this GitHub repo, including it in your project (say the include folder) and then adding the following lines to the CMakeLists.txt file,
```sh
#some code
include_directories("${CMAKE_SOURCE_DIR}/include/MorfiConvert")
add_subdirectory(MorfiConvert)
#some more code
set(..
add_executable(…
target_link_libraries(Your_Project_Name morficonvert)
```

Here “Your_Project_Name” is the name of the project that is to use MorfiConvert. 

Afterwards include the MorfiConvert and GeneData (A class within MorfiConvert as an intermediary data state) header files.
```sh
#include "include/MorfiConvert/src/morficonvert.h"
#include "include/MorfiConvert/src/GeneData.h"
...
```

or

MorfiConvert is in cppan as well.
Make a "cppan.yml" file with,
```sh
dependencies:
    pvt.damitha_lenadora.morficonvert: 1.0.2
```
Then go to the directory your project is in and build it by (you'll need cppan downloaded and in your path variables),,
```sh
cppan
mkdir build && cd build
cmake ..
cmake --build . --config Release
```
or

If your going with CMake to add it to your project, just,
```sh
cppan
```
And update the CMakeLists.txt file.
```sh
cmake_minimum_required(VERSION 3.11)
project(Someproject)

set(CMAKE_CXX_STANDARD 17)

# Use solution folders.
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMake Targets")

# Output directory settings
set(output_dir ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${output_dir})

if (MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

add_subdirectory(.cppan)

#Add whatever you want to the project
...

target_link_libraries(Someproject
        pvt.damitha_lenadora.morficonvert
        )
```
To use it like,
```sh
#include <morficonvert.h>
#include <GeneData.h>
```


# Functions

## Conversion Functions
### Text file to GeneData Conversion
```sh
std::vector<GeneData> getInData(std::string path, std::string origin, int thrdCnt = 10)
```
Converts an input data file into an intermediary format. 
Here path specifies the file path to the input data file, origin specifies the input file format and thrdCnt gives the user the ability to specify the number of threads that the function should run on (default 10). Supported formats are,
- NCBI-BLAST – .txt, .xml, .json
- BLAST – .txt
- MSA – .txt
- FASTA – Any so long as the data in it is FASTA
- GCG – Any so long as the data in it is GCG
- EMBL – Any so long as the data in it is EMBL
- GenBank – Any so long as the data in it is GenBank
- PIR – Any so long as the data in it is PIR
- NBRF – Any so long as the data in it is NBRF
- PHYLIP – Any so long as the data in it is PHYLIP
- Clustal – .txt
- MSF – .txt
- FASTA-Report – .txt
- FIND - MorfiConvert tries to figure out which type of file input format from the above is given and then converts

Ex: - 
```sh
std::vector<GeneData> inpData = morfiConvert::getInData(“F:\\SomeFolder\\Test1.txt”, “NCBI-BLAST”); 
```
### Image file to GeneData Conversion
```sh
std::vector<GeneData> getInImageData(std::vector<std::string> genes, std::string path)
```

Converts an input image file in the MSA format into GeneData objects.
Here, genes specify the descriptions of the said gene sequences in the image file. Supported formats are,

- .bmp

The image you have is not .bmp? Try converting it. You can find loads of sites that do image to .bmp conversions.

Ex: - 
```sh
std::vector<GeneData> inpData = morfiConvert::getImageData(geneDescrips, “F:\\SomeFolder\\img.bmp”); 
```
`But keep in mind that the data in the image needs to be trained first. This is done via tainImageData() passing the gene sequence as a vector set of strings and the image file path as parameters.`
```sh
std::vector<GeneData> trainImageData(std::vector<std::string> genes, std::string path)
```
### GeneData to the required output format (as string objects)
```sh
std::vector<std::string> giveOutData(std::vector<GeneData> res, std::string outType, int thrdCnt = 10);
```

Converts GeneData to a set strings of the requested format. Such requested formats are,
- FASTA – FASTA object set with the description in the first line followed by the gene
- NCBI-IDs – The NCBI identifiers
- GCG
- EMBL
- GenBank
- PIR
- NBRF
- PHYLIP

Ex:-
```sh
std::vector<std::string> outData = morfiConvert::giveOutData(res,”FASTA”);
```

###  GeneData to the required output format (as a file)
```sh
void writeOutData(std::vector<GeneData> res, std::string path, std::string outType, int thrdCnt = 10);
```
Converts GeneData to a file of the requested format. Such requested formats are,
- FASTA
- GCG
- GenBank
- PIR
- NBRF
- PHYLIP

Ex:-
```sh
morfiConvert::writeOutData(res,”FASTA”);
```

## Intermediary Functions
### Identify the given data file format
```sh
std::string identify(std::string path);
```
Identify the data format of the input file. Forms identified are,
- NCBI-BLAST
- GCG
- FASTA
- EMBL 
- GenBank
- PIR / NBRF
- Clustal
- FASTA-Report

Ex:-
```sh
std::string inpType = morfiConvert:: identify (“F:\\SomeFolder\\Test1.txt”); 
```

### Create a local copy of FASTA database files
```sh
void configLocal(std::string path, std::string outPath = "");
```
Creates a local data file from NCBI FASTA database files. Primarily done to give MorfiConvert the capability to complete partial gene sequences. An output location can also be provided via outPath. If not given the local files would be created in the same folder as the C++ executable. 
Ex:-
```sh
morfiConvert:: configLocal (“F:\\SomeFolder\\Test1.txt”);
```

### Get the full gene sequence (if in the local files)
```sh
std::string fullGene(std::string ncbiID, std::string dbName = "all");
```
Search the local data files and find the matching gene via the NCBI identifier. The db parameter refers to the name of the database file specified to find a particular record.
>Ex:- If the FASTA database file was named “env_nr” used in configLocal function, and the records of this said file needed to be searched the db parameter should be “env_nr”.  

A point to note would be that if the db parameter was specified as “all” then all available local data files would be searched and that if no matching record was found an empty string would be returned.
Ex:-
```sh
morfiConvert::completeGene(“ACR012.1”,”all”);
```

### Get the full sequences of a GeneData set (if in the local files)
```sh
std::vector<GeneData> completeGenes(std::vector<GeneData> geneSet);
```
The multi object format of completeGene. The differences to note in this function to that of completeGene would be that if a matching record was not found or if the record in question was found but the gene sequence was smaller, then the current gene would not be changed, as well as having all local database files searched.
Ex:-
```sh
morfiConvert::competeGenes(geneSet);
```

### Gene fragment via index
```sh
GeneData indexFragment(GeneData inpData, int strtIndex, int endIndex);
```
Take only the relevant part of the gene sequence. Done by the location of the gene from the start (starting form 0) and the ending location.
Ex:-
```sh
morfiConvert::indexFragment(gene, 0, 5);
```

### Gene fragment via string
```sh
GeneData stringFragment(GeneData inpData, std::string strtString, std::string endString);
```
Take only the relevant part of the gene sequence. Gives the section of the gene from the first matching string found in the current gene to the last matching ending string. 
Ex:-
```sh
morfiConvert::indexFragment(gene, “MSA”, “YPP”);
```

### Gene set fragment via index
```sh
std::vector<GeneData> multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex);
```
The multi object format of indexFragmet.

### Gene set fragment via string
```sh
std::vector<GeneData> multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString);
```
The multi object format of stringFragment.

## Complete Conversions 
### Convert from file to final form as objects
```sh
std::vector<std::string> objConvert(std::string inPath, std::string origin, std::string outType, int thrdCnt=10);
```
The complete form of file to object conversion where the output is produced directly from the input.

### Convert from file to final form a file
```sh
void fileConvert(std::string inPath, std::string outPath, std::string origin, std::string outType, int thrdCnt=10);
```
The complete form of file to file conversion where the output is produced directly from the input.

## Testing
MorfiConvert has tests created using gtest (in the github repo, not in the cppan one). Feel free to try them out. `However be warned that if the test directory is the same as the build directory then locally saved data can be erased. This is due to MorfiTest cleaning up after itself.` 

License
----

MIT