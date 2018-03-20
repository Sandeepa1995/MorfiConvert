# MorfiConvert : Bioinformatics Data Format Converter

MorfiConvert is a bioinformatics data format converter created to specialize in converting the outputs of one bioinformatics process to a standard format that can be used as an input for another bioinformatics process.

# Using MorfiConvert 
MorfiConvert is uses CMake to build its file to give itself platform independency. Thus to use MorfiConvert you will need to,
1.	Have a C++ compiler and install CMake
2.	Build the project via entering the following in the command-line within the MorfiConvert folder,
```sh
mkdir build && cd build
cmake ..
```
or 
Using the MorfiConvert in a project where CMake is used via adding the following lines to the CMake file,
```sh
#some code
include_directories("${CMAKE_SOURCE_DIR}/MorfiConvert/src")
link_directories(MorfiConvert/src)
add_subdirectory(MorfiConvert/src)
#some more code
set(..
add_executable(…
target_link_libraries(Your_Project_Name morficonvert)
```
*Here “Your_Project_Name” is the name of the project that is to use MorfiConvert.

# Functions
## Conversion Functions
### getInData
```sh
std::vector<GeneData> getInData(std::string path, std::string origin)
```
Converts an input data file into an intermediary format (A class within MrofiConvert – GeneData Class).
Here path specifies the file path to the input data file and origin specifies the input file format. Supported formats are,
- NCBI-BLAST – .txt, .xml, .json
- BLAST – .txt
- MSA – .txt
- FASTA – Any so long as the data in it is FASTA
- FIND – .txt

By giving the origin as FIND MorfiConvert tries to figure out which type of file input format is given in the file mentioned. 
Ex: - 
```sh
std::vector<GeneData> inpData = morfiConvert::getInData(“F:\\SomeFolder\\Test1.txt”, “NCBI-BLAST”); 
```

### giveOutData
```sh
std::vector<std::string> giveOutData(std::vector<GeneData> res, std::string outType);
```

Converts GeneData to a set strings of the requested format. Such requested formats are,
- FASTA – FASTA object set with the description in the first line followed by the gene
- NCBI-IDs – The NCBI identifiers

Ex:-
```sh
std::vector<std::string> outData = morfiConvert::giveOutData(res,”FASTA”);
```

writeOutData
void writeOutData(std::vector<GeneData> res, std::string path, std::string outType);
Converts GeneData to a file of the requested format. Such requested formats are,
FASTA
ex:-
morfiConvert::writeOutData(res,”FASTA”);

## Intermediary Functions
### identify
```sh
std::string identify(std::string path);
```
Identify the data format of the input file. Specifically, if the format is NCBI-BLAST.
Ex:-
```sh
std::string inpType = morfiConvert:: identify (“F:\\SomeFolder\\Test1.txt”); 
```

### configLocal
```sh
void configLocal(std::string path);
```
Create a local data file from NCBI FASTA db files. Primarily done to give MorfiConvert the capability to complete partial gene sequences. 
Ex:-
```sh
morfiConvert:: configLocal (“F:\\SomeFolder\\Test1.txt”);
```

### completeGene
```sh
std::string completeGene(std::string ncbiID, std::string db);
```
Search the local data files and find the matching gene via the NCBI identifier. The “db” parameter refers to the name of the database file specified to find a particular record. 
>Ex:- If the FASTA database file was named “env_nr” used in configLocal function, and the records of this said file needed to be searched the db parameter should be “env_nr”. 

A point to note would be that if the db parameter was specified as “all” then all available local data files would be searched and that if no matching record was found an empty string would be returned.
Ex:-
```sh
morfiConvert::completeGene(“ACR012.1”,”all”);
```

### completeGenes
```sh
std::vector<GeneData> completeGenes(std::vector<GeneData> geneSet);
```
The multi object format of completeGene. Differences to note in this function to that of completeGene would be that if a matching record was not found or if the record in question was found but the gene was smaller, then the current gene would not be changed and that all local database files would be searched.
Ex:-
```sh
morfiConvert::competeGenes(geneSet);
```

### indexFragmet
```sh
GeneData indexFragment(GeneData inpData, int strtIndex, int endIndex);
```
Take only the relevant part of the gene sequence. Done by the location of the gene from the start and the length of the final sequence. 
Ex:-
```sh
morfiConvert::indexFragment(gene, 0, 5);
```

### stringFragment
```sh
GeneData stringFragment(GeneData inpData, std::string strtString, std::string endString);
```
Take only the relevant part of the gene sequence. Done by the first matching string found in the current gene to the last matching ending gene portion of the final sequence. 
Ex:-
```sh
morfiConvert::indexFragment(gene, “MSA”, “YPP”);
```

### multiIndexFragment
```sh
std::vector<GeneData> multiIndexFragment(std::vector<GeneData> inpData, int strtIndex, int endIndex);
```
The multi object format of indexFragmet.

### multiStringFragment
```sh
std::vector<GeneData> multiStringFragment(std::vector<GeneData> inpData, std::string strtString, std::string endString);
```
The multi object format of stringFragment.

## Complete Conversions 
### objConvert
```sh
std::vector<std::string> objConvert(std::string inPath, std::string origin, std::string outType);
```
The complete form of file to object conversion where the output is produced directly from the input.

### fileConvert
```sh
void fileConvert(std::string inPath, std::string outPath, std::string origin, std::string outType);
```
The complete form of file to file conversion where the output is produced directly from the input.

License
----

MIT