//
// Created by Damitha on 3/21/2018.
//
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <exception>
#include <omp.h>

#include "morficonvert.h"
#include "supplimentMorfi.h"
#include "compressGene.h"

#include "gtest/gtest.h"

//**********************************************************************************************************************
//Predefined paths
//**********************************************************************************************************************

//!!!!!!!!!!!!!!!!!!!!Path will need to be changed - Some files are in the Tests subdirectory!!!!!!!!!!!!!!!!!!!!!!!!!!!

//Path where existing library is built
std::string prePath = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\untitled6\\cmake-build-debug\\bin\\";

//Path to the FASTA data file from NCBI
std::string datapath = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Data\\env_nr\\env_nr";

//Incorrect path
std::string wrongpath = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Data\\env_nr\\env.txt";

//Any existing file other than a FASTA data file
std::string wrongfile = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\1 Template for Project Proposal.doc";

//NCBI .txt file for complete valid data
std::string pathNCBITXTValid = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Valid Total (NCBI)\\Test1.txt";

//NCBI .xml file for complete valid data
std::string pathNCBIXMLValid = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Valid Total (NCBI)\\Test1.xml";

//NCBI .json file for complete valid data
std::string pathNCBIJSONValid = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Valid Total (NCBI)\\Test1.json";

//Custom MSA
std::string pathCusMSA = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\MSA Text\\Test1.txt";

//Custom BLAST
std::string pathCusBLAST = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Valid (Custom)\\Test1.txt";

//FASTA
std::string pathFASTA = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\FASTA\\Test1.txt";

//NCBI .txt file for partially valid data
std::string pathNCBITXTInvalid = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Invalid\\Test1.txt";

//NCBI .xml file for partially valid data
std::string pathNCBIXMLInvalid = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Invalid\\Test1.xml";

//NCBI .json file for partially valid data
std::string pathNCBIJSONInvalid = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Invalid\\Test1.json";

//GCG
std::string pathGCG = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\EBI Data\\GCG.txt";

//EMBL
std::string pathEMBL = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\EBI Data\\EMBL.txt";

//GenBank
std::string pathGenBank = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\EBI Data\\GenBank.txt";

//MSF
std::string pathMSF = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\EBI Data\\MSF.txt";

//PHYLIP
std::string pathPHYLIP = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\EBI Data\\PHYLIPD.txt";

//PIR
std::string pathPIR = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\EBI Data\\PIR.txt";

//Clustral
std::string pathClus = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\EBI Data\\Clustral.txt";

//FASRA Report
std::string pathFR = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\EBI Data\\FASTARes.txt";

//MSA Image
std::string msaImage = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\download.bmp";

//**********************************************************************************************************************
//Data Test
//**********************************************************************************************************************

//Check local gene storage.
bool verifyLocalData(std::string s) {
    int i = 0;
    while (s[i] != '.' && s[i] != '_') {
        if (!(isalnum(s[i]))) {
            return false;
        }
        i += 1;
    }
    i += 1;
    while (s[i] != ' ') {
        if (!(isalnum(s[i]))) {
            return false;
        }
        i += 1;
    }
    std::string gene = CompressGene::decodeGene(s.substr(i + 1));
    int j = 0;
    while (j < gene.length()) {
        if (!(isalpha(gene[j]))) {
            return false;
        }
        j += 1;
    }
    return true;
}

//Verify the data locally stored
TEST(data_test, test_verify_local_genes) {
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(prePath + "points.morfi");

    std::string s;
    std::string sGene;

    if (infile.good()) {
        while (getline(infile, s)) {
            //String without the final pointer value
            std::string shaveLastPoint = s.substr(0, s.find_last_of(" "));

            //Get the file path of the local database file
            std::string dbPath = shaveLastPoint.substr(0, shaveLastPoint.find_last_of(" "));

            std::ifstream datafile;
            datafile.open(dbPath);

            while (getline(datafile, sGene)) {
                EXPECT_TRUE(verifyLocalData(sGene));
            }
        }
    } else {
        EXPECT_TRUE(false);
    }

    std::cout << "Data test - test_verify_local_genes : Complete" << std::endl;
}

//Check local character recognition
bool verifyLocalChar(std::string s) {
    int i = 0;
    while (s[i] != ' ') {
        if (!(isdigit(s[i]))) {
            return false;
        }
        i += 1;
    }
    i += 1;
    while (i < s.length()) {
        if (s[i] != '0' && s[i] != '1') {
            return false;
        }
        i += 1;
    }
    return true;
}

//Verify the locally saved character recognition data
TEST(data_test, test_verify_local_char_rec) {
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open(prePath + "Character patterns.morfi");

    std::string s;

    if (infile.good()) {
        while (getline(infile, s)) {
            EXPECT_TRUE(verifyLocalChar(s));
        }
    } else {
        EXPECT_TRUE(false);
    }

    std::cout << "Data test - test_verify_local_char_rec : Complete" << std::endl;
}

//**********************************************************************************************************************
//Function Test
//**********************************************************************************************************************

//Test local storage (No errors)
TEST(function_test, test_local_storage_correct) {
    try {
        morfiConvert::configLocal(datapath);

        //Input file stream
        std::ifstream infile;

        //Open file
        infile.open("points.morfi");

        std::string s;
        std::string sGene;

        std::string geneFound = "";

        if (infile.good()) {
            while (getline(infile, s)) {
                //String without the final pointer value
                std::string shaveLastPoint = s.substr(0, s.find_last_of(" "));

                //Get the file path of the local database file
                std::string dbPath = shaveLastPoint.substr(0, shaveLastPoint.find_last_of(" "));

                std::ifstream datafile;
                datafile.open(dbPath);

                while (getline(datafile, sGene)) {
                    EXPECT_EQ(true, verifyLocalData(sGene));
                }
                datafile.close();
                remove(dbPath.c_str());
            }
            remove("points.morfi");
        } else {
            EXPECT_EQ(true, false);
        }
    }
    catch (std::exception &e) {
        EXPECT_EQ("Expect exception not thrown", "Exception thrown.");
    }

    std::cout << "Function Test - test_local_storage_correct : Complete" << std::endl;
}

//Test local storage (Wrong path)
TEST(function_test, test_local_storage_errpath) {
    try {
        morfiConvert::configLocal(wrongpath);
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_local_storage_errpath : Complete" << std::endl;
}

//Test local storage (Wrong file)
TEST(function_test, test_local_storage_wrongfile) {
    try {
        morfiConvert::configLocal(wrongfile);
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in the file type given. Did not contain data in the FASTA format.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_local_storage_wrongfile : Complete" << std::endl;
}

TEST(function_test, test_gene_complete) {
    morfiConvert::configLocal(datapath);

    std::string expResult = "MTGQRIGYIRVSTFDQNPERQLEGVKVDRAFSDKASGKDVKRPQLEALISFARTGDTVVVHSMDRLARNLDDLRRIVQTLTQRGVHIEFVKEHLSFTGEDSPMANLMLSVMGAFAEFERALIRERQREGIALAKQRGAYRGRKKSLSSERIAELRQRVEAGEQKTKLAREFGISRETLYQYLRTDQ";
    std::string testval1 = morfiConvert::fullGene("EBT52160.1", "env_nr");
    EXPECT_EQ(testval1, expResult);
    std::string testval2 = morfiConvert::fullGene("EBT52160.1", "all");
    EXPECT_EQ(testval2, expResult);
    std::string testval3 = morfiConvert::fullGene("EBT52160.1", "env_np");
    EXPECT_EQ(testval3, "");
    std::string testval4 = morfiConvert::fullGene("EBT521650.1", "env_nr");
    EXPECT_EQ(testval4, "");

    std::string expResult2 = "AGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMSXYQTPLFVWSVLITAVLLLLSLPVLATGITMLLXDRNLNTTFFXPXXGGDP";
    std::string testval5 = morfiConvert::fullGene("ARO40711.1", "env_nr");
    EXPECT_EQ(testval5, expResult2);

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open("points.morfi");

    std::string s;
    std::string sGene;

    std::string geneFound = "";

    if (infile.good()) {
        while (getline(infile, s)) {
            //String without the final pointer value
            std::string shaveLastPoint = s.substr(0, s.find_last_of(" "));

            //Get the file path of the local database file
            std::string dbPath = shaveLastPoint.substr(0, shaveLastPoint.find_last_of(" "));

            remove(dbPath.c_str());
        }
        remove("points.morfi");
    }

    std::cout << "Function Test - test_gene_complete : Complete" << std::endl;

}

std::vector<GeneData> getEBIDataTest() {
    std::vector<GeneData> outRes;
    GeneData r1("BAA20512.1",
                "ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGGGCTAAGATCAGCCCCAAAGCCGATGATATCGGCGCTGAAGCTCTCGGCAGAATGCTGACCGTCTACCCTCAGACCAAGACCTACTTCGCTCACTGGGATGACCTGAGCCCTGGGTCCGGTCCTGTGAAGAAGCATGGCAAGGTTATCATGGGTGCAGTGGCCGATGCCGTTTCAAAAATAGACGACCTTGTGGGAGGTCTGGCCTCCCTGAGCGAACTTCATGCTTCCAAGCTGCGTGTTGACCCGGCCAACTTCAAGATCCTCGCACACAATGTCATCGTGGTCATCGGCATGCTCTTCCCTGGAGACTTCCCCCCAGAGGTTCACATGTCAGTTGACAAGTTTTTCCAGAACTTGGCTCTGGCTCTCTCTGAGAAGTACCGCTAA");
    outRes.push_back(r1);
    GeneData r2("CAA23748.1",
                "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAAGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCTTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAA");
    outRes.push_back(r2);
    GeneData r3("CAA24095.1",
                "ATGGTGCTCTCTGGGGAAGACAAAAGCAACATCAAGGCTGCCTGGGGGAAGATTGGTGGCCATGGTGCTGAATATGGAGCTGAAGCCCTGGAAAGGATGTTTGCTAGCTTCCCCACCACCAAGACCTACTTTCCTCACTTTGATGTAAGCCACGGCTCTGCCCAGGTCAAGGGTCACGGCAAGAAGGTCGCCGATGCGCTGGCCAGTGCTGCAGGCCACCTCGATGACCTGCCCGGTGCCTTGTCTGCTCTGAGCGACCTGCATGCCCACAAGCTGCGTGTGGATCCCGTCAACTTCAAGCTCCTGAGCCACTGCCTGCTGGTGACCTTGGCTAGCCACCACCCTGCCGATTTCACCCCCGCGGTACATGCCTCTCTGGACAAATTCCTTGCCTCTGTGAGCACCGTGCTGACCTCCAAGTACCGTTAA");
    outRes.push_back(r3);
    GeneData r4("CAA28435.1",
                "ATGTCTCTGACCAGGACTGAGAGGACCATCATCCTGTCCCTGTGGAGCAAGATCTCCACACAGGCAGACGTCATTGGCACCGAGACCCTGGAGAGGCTCTTCTCCTGCTACCCGCAGGCCAAGACCTACTTCCCGCACTTCGACCTGCACTCGGGCTCCGCGCAGCTGCGCGCGCACGGCTCCAAGGTGGTGGCCGCCGTGGGCGACGCGGTCAAGAGCATCGACAACGTGACGAGCGCGCTGTCCAAGCTGAGCGAGCTGCACGCCTACGTGCTGCGCGTGGACCCGGTCAACTTCAAGTTCCTGTCCCACTGCCTGCTGGTCACGTTGGCCTCGCACTTCCCCGCCGACTTCACGGCCGACGCGCACGCCGCCTGGGACAAGTTCCTGTCCATCGTGTCCGGCGTCCTGACGGAGAAGTACCGCTGA");
    outRes.push_back(r4);

    return outRes;
}

std::vector<GeneData> getEBIDataFASTARepTest() {
    std::vector<GeneData> outRes;
    GeneData r1("EX881837.1",
                "AATCTAACTTGAGAAAAAGAAGACGCAGCAATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGGGCTAAGATCAGCCCCAAAGCCGATGATATCGGCGCTGAAGCTCTCGGCAGAATGCTGACCGTCTACCCTCAGACCAAGACCTACTTCGCTCACTGGGCTGACCTGAGCCCTGGGTCCGGTCCTGTGAAGAAGCATGGCAAGGTTATCATGGGTGCAGTCGGCGATGCCGTTTCAAAAATAGACGACCTTGTGGGAGGTCTGGCCTCCCTGAGCGAACTTCATGCTTCCAAGCTGCGTGTTGACCCGGCCAACTTCAAGATCCTCGCACACAATGTCATCGTGGTCATCGGCATGCTCTTCCCTGGAGACTTCCCCCCAGAGGTTCACATGTCAGTTGACAAGTTTTTCCAGAACTTGGCTCTGGCTCTCTCTGAGAAGTACCGCTAAATTTCCGGTGGGCATTCATACGCGACACTGTGCGGCACTCCTGACCAACTCAAGTGTGATGTCTGAATAAAATTTCTC");
    outRes.push_back(r1);
    GeneData r2("EX882361.1",
                "AATCTAACTTGAGAAAAAGAAGACGCGGCAATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGGGCTAAGATCAGCCCCAAAGCCGATGATATCGGCGCTGAAGCCCTCGGCAGAATGCTGACCGTCTACCCTCAGACCAAGACCTACTTCGCTCACTGGGCTGACCTGAGCCCTGGGTCCGGTCCTGTGAAGAAGCATGGCAAGGTTATCATGGGTGCAGTCGGCGATGCCGTTTCAAAAATAGACGACCTTGTGGGAGGTCTGGCTGCCCTGAGCGAACTTCATGCTTTCAAGCTGCGTGTTGACCCGGCCAACTTCAAGATCCTCGCACACAATGTCATCGTGGTCATCGGCATGCTCTACCCTGGAGACTTCCCCCCAGAGGTTCACATGTCAGTTGACAAGTTTTTCCAGAACTTGGCTCTGGCTCTCTCTGAGAAGTACCGCTAAATTTCCGGTGGGCATTCATACGCGACACTGTGCGGCACTCCTGACCAACTCAAGTGTGATGTCTGAATAAAATTTCTC");
    outRes.push_back(r2);

    return outRes;
}

//Test identifiyer extraction
TEST(function_test, test_get_identifier) {
    EXPECT_EQ(SupMorfi::getIdentifier("BAA20512.1"), "BAA20512.1");
    EXPECT_EQ(SupMorfi::getIdentifier("ENA|BAA20512|BAA20512.1 Cyprinus carpio (common carp) alpha-globin"),
              "BAA20512.1");
    EXPECT_EQ(SupMorfi::getIdentifier("BAA20512.1 Cyprinus carpio (common carp) alpha-globin"), "BAA20512.1");

    std::cout << "Function Test - test_get_identifier : Complete" << std::endl;
}

std::vector<GeneData> getGeneDataNCBIBLASTTest() {
    std::vector<GeneData> outRes;
    GeneData r1("AAD44166.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVILGLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGXIENY");
    outRes.push_back(r1);
    GeneData r2("BAA25010.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGMIENY");
    outRes.push_back(r2);
    GeneData r3("YP_626379.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSNSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGMIENY");
    outRes.push_back(r3);
    GeneData r4("AAD44169.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEHPYIIIGQMASILYFSIILAFLPIAGMIENY");
    outRes.push_back(r4);
    GeneData r5("BAA25017.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEHPYIIIGQMASILYFSIILAFLPIAGMIENY");
    outRes.push_back(r5);

    return outRes;
}

//Test NCBI BLAST data extraction in .txt file
TEST(function_test, test_file_2_obj_NCBI_txt_valid) {
    //NCBI BLAST data with no errors in the file
    std::vector<GeneData> compareSet = getGeneDataNCBIBLASTTest();

    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");
    int geneMatches = 0;
    int desMatches = 0;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < test1.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(test1.at(j).getDetails()).c_str(),
                       compareSet.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(test1.at(j).getGene().c_str(), compareSet.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }

    std::cout << "Function Test - test_file_2_obj - BLAST .txt all valid : Complete" << std::endl;

    EXPECT_EQ(desMatches, 5);
    EXPECT_EQ(geneMatches, 5);
    EXPECT_EQ(test1.size(), 100);

}

//Test NCBI BLAST data extraction in .xml file
TEST(function_test, test_file_2_obj_NCBI_xml_valid) {
    //NCBI BLAST data with no errors in the file
    std::vector<GeneData> compareSet = getGeneDataNCBIBLASTTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> test2 = morfiConvert::getInData(pathNCBIXMLValid, "NCBI-BLAST");
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < test2.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(test2.at(j).getDetails()).c_str(),
                       compareSet.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(test2.at(j).getGene().c_str(), compareSet.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 5);
    EXPECT_EQ(geneMatches, 5);
    EXPECT_EQ(test2.size(), 100);

    std::cout << "Function Test - test_file_2_obj - BLAST .xml all valid: Complete" << std::endl;
}

//Test NCBI BLAST data extraction in .json file
TEST(function_test, test_file_2_obj_NCBI_json_valid) {
    //NCBI BLAST data with no errors in the file
    std::vector<GeneData> compareSet = getGeneDataNCBIBLASTTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> test3 = morfiConvert::getInData(pathNCBIJSONValid, "NCBI-BLAST");
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < test3.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(test3.at(j).getDetails()).c_str(),
                       compareSet.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(test3.at(j).getGene().c_str(), compareSet.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 5);
    EXPECT_EQ(geneMatches, 5);
    EXPECT_EQ(test3.size(), 100);

    std::cout << "Function Test - test_file_2_obj - BLAST .json all valid: Complete" << std::endl;


}

std::vector<GeneData> getGeneDataMSATextTest() {
    std::vector<GeneData> outRes;
    GeneData r1("AGS56207.1",
                "MSAVQSAITTAEQQHRIYKAFDWVVIPIVILVITGAFHIHAMLTMGDWDFWIDWKDRRWWVTVTPIALIFFPAALHYLLWTNFRLPFGATLSCSVLLIGEWANRTANFVGWA");
    outRes.push_back(r1);
    GeneData r2("AGS56208.1",
                "MSASQSAVRSRAEAVAVSRTFDWMILFVLFTAVLGGYHIHFMLTGGDWDFWSDWKDRRLWVTVVPIVAITFPAAVQACLWWRLPIGATLCILALLLGEWINRYINFWGWTYF");
    outRes.push_back(r2);
    GeneData r3("AGS56209.1",
                "MGVVRSAITTAEQQHRVYKAFDWIVIPIVVLVIIGAFHIHFMLTAGDWDFWIDWKDRRWWVTITPITLIFFPAALHYLFWQRLPFGATICCAGLLAGEWVSRVVNFVGWAYY");
    outRes.push_back(r3);
    GeneData r4("AGS56210.1",
                "MSLTAEEAKLLKAISVPNSDAAKIVRNFDWIVVVVLLFAVTAAFHIHFMLTAGDWDFWVDWKDRMYWVTLTPMILVIIPAAYPIMHKLGLPIAGTVCGLCLVIGEWVSRISG");
    outRes.push_back(r4);
    GeneData r5("AGS56211.1",
                "MAYSDAITQKLNMSTNEERTYRLFDLIIIPIVVVVLGAAFHLHYMLTAGDWDFWVDWKDRLYWPLISSVSFVIFAGAVQGLWRTIRMPIGMTVAAIGLIAMEWISRYVNFYE");
    outRes.push_back(r5);
    GeneData r6("AGS56212.1",
                "MSSTQSAVRSHAEAVQVSRTIDYLGLFILFFVILGGFHVHAMLTMGDWDFWSDWKDRRLWVTVTPIMLVTFPAAVQAIVWEIGFGATLCCISLVLGEWINRYFNFWGWTYFP");
    outRes.push_back(r6);

    return outRes;
}

//Test MSA custom data extraction
TEST(function_test, test_file_2_obj_MSA_Custom) {
    //Custom MSA data
    std::vector<GeneData> compareSetMSA = getGeneDataMSATextTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> test4 = morfiConvert::getInData(pathCusMSA, "MSA");
    for (int i = 0; i < 6; i++) {
        EXPECT_EQ(test4.at(i).getGene(), compareSetMSA.at(i).getGene());
        EXPECT_EQ(test4.at(i).getDetails(), compareSetMSA.at(i).getDetails());
    }
    EXPECT_EQ(test4.size(), 6);

    std::cout << "Function Test - test_file_2_obj - Custom MSA : Complete" << std::endl;

}

std::vector<GeneData> getGeneDataCustomBLASTTest() {
    std::vector<GeneData> outRes;
    GeneData r1("AAD44166.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVILGLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGXIENY");
    outRes.push_back(r1);
    GeneData r2("ABD44166.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVILGLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGXIENY");
    outRes.push_back(r2);
    GeneData r3("ACD44166.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVILGLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGXIENY");
    outRes.push_back(r3);
    GeneData r4("ATX68788.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTDLVEWIWGGFSVDKATLNRFFALHFILPFTMIALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILFLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALLLSILILGIMPLLHTSKHRSMMLRPLSQVLFWTLATDLLMLTWIGSQPVEYPYIIIGQMASILYFSIILVFLPIAGMIENY");
    outRes.push_back(r4);
    GeneData r5("ATY68788.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTDLVEWIWGGFSVDKATLNRFFALHFILPFTMIALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILFLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALLLSILILGIMPLLHTSKHRSMMLRPLSQVLFWTLATDLLMLTWIGSQPVEYPYIIIGQMASILYFSIILVFLPIAGMIENY");
    outRes.push_back(r5);

    return outRes;
}

//Test BLAST custom data extraction
TEST(function_test, test_file_2_obj_BLAST_Custom) {
    //Custom BLAST data
    std::vector<GeneData> compareSetCusBLAST = getGeneDataCustomBLASTTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> test5 = morfiConvert::getInData(pathCusBLAST, "BLAST");
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < test5.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(test5.at(j).getDetails()).c_str(),
                       compareSetCusBLAST.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(test5.at(j).getGene().c_str(), compareSetCusBLAST.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 5);
    EXPECT_EQ(geneMatches, 5);
    EXPECT_EQ(test5.size(), 5);

    std::cout << "Function Test - test_file_2_obj - Custom BLAST : Complete" << std::endl;

}

//Test FASTA data extraction
TEST(function_test, test_file_2_obj_FASTA) {
    //NCBI BLAST data with no errors in the file
    std::vector<GeneData> compareSet = getGeneDataNCBIBLASTTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> test6 = morfiConvert::getInData(pathFASTA, "FASTA");
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < test6.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(test6.at(j).getDetails()).c_str(),
                       compareSet.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(test6.at(j).getGene().c_str(), compareSet.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 5);
    EXPECT_EQ(geneMatches, 5);
    EXPECT_EQ(test6.size(), 100);

    std::cout << "Function Test - test_file_2_obj - FASTA : Complete" << std::endl;

}

//Test type identified data extraction
TEST(function_test, test_file_2_obj_find) {
    //NCBI BLAST data with no errors in the file
    std::vector<GeneData> compareSet = getGeneDataNCBIBLASTTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> test7 = morfiConvert::getInData(pathNCBITXTValid, "FIND");
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < test7.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(test7.at(j).getDetails()).c_str(),
                       compareSet.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(test7.at(j).getGene().c_str(), compareSet.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 5);
    EXPECT_EQ(geneMatches, 5);
    EXPECT_EQ(test7.size(), 100);

    std::cout << "Function Test - test_file_2_obj - Find and extract - BLAST .txt : Complete" << std::endl;

}

//Test type GCG data extraction
TEST(function_test, test_file_2_obj_GCG) {
    //Data for the GCG, EMBL, GenBank, MSF, PHYLIP, PIR and Clustral data formats taken from the EBI website
    std::vector<GeneData> compareSetEBI = getEBIDataTest();

    int geneMatches = 0;
    int desMatches = 0;

    //EBI data
    std::vector<GeneData> testGCG = morfiConvert::getInData(pathGCG, "GCG");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < testGCG.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(testGCG.at(j).getDetails()).c_str(),
                       compareSetEBI.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(testGCG.at(j).getGene().c_str(), compareSetEBI.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 4);
    EXPECT_EQ(geneMatches, 4);
    EXPECT_EQ(testGCG.size(), 4);

    std::cout << "Function Test - test_file_2_obj - GCG : Complete" << std::endl;
}

//Test type EMBL data extraction
TEST(function_test, test_file_2_obj_EBML) {
    //Data for the GCG, EMBL, GenBank, MSF, PHYLIP, PIR and Clustral data formats taken from the EBI website
    std::vector<GeneData> compareSetEBI = getEBIDataTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> testEMBL = morfiConvert::getInData(pathEMBL, "EMBL");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < testEMBL.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(testEMBL.at(j).getDetails()).c_str(),
                       compareSetEBI.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(testEMBL.at(j).getGene().c_str(), compareSetEBI.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 4);
    EXPECT_EQ(geneMatches, 4);
    EXPECT_EQ(testEMBL.size(), 5);

    std::cout << "Function Test - test_file_2_obj - EMBL : Complete" << std::endl;
}

//Test type GenBank data extraction
TEST(function_test, test_file_2_obj_Genbank) {
    //Data for the GCG, EMBL, GenBank, MSF, PHYLIP, PIR and Clustral data formats taken from the EBI website
    std::vector<GeneData> compareSetEBI = getEBIDataTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> testGenBank = morfiConvert::getInData(pathGenBank, "GenBank");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < testGenBank.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(testGenBank.at(j).getDetails()).c_str(),
                       compareSetEBI.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(testGenBank.at(j).getGene().c_str(), compareSetEBI.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 4);
    EXPECT_EQ(geneMatches, 4);
    EXPECT_EQ(testGenBank.size(), 6);

    std::cout << "Function Test - test_file_2_obj - GenBank : Complete" << std::endl;
}

//Test type MSF data extraction
TEST(function_test, test_file_2_obj_MSF) {
    //Data for the GCG, EMBL, GenBank, MSF, PHYLIP, PIR and Clustral data formats taken from the EBI website
    std::vector<GeneData> compareSetEBI = getEBIDataTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> testMSF = morfiConvert::getInData(pathMSF, "MSF");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < testMSF.size(); j++) {

            if (strcmp(SupMorfi::getIdentifier(testMSF.at(j).getDetails()).c_str(),
                       compareSetEBI.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(testMSF.at(j).getGene().c_str(), compareSetEBI.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 4);
    EXPECT_EQ(geneMatches, 4);
    EXPECT_EQ(testMSF.size(), 4);

    std::cout << "Function Test - test_file_2_obj - MSF : Complete" << std::endl;
}

//Test type PHYLIP data extraction
TEST(function_test, test_file_2_obj_PHYLIP) {
    //Data for the GCG, EMBL, GenBank, MSF, PHYLIP, PIR and Clustral data formats taken from the EBI website
    std::vector<GeneData> compareSetEBI = getEBIDataTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> testPHYLIP = morfiConvert::getInData(pathPHYLIP, "PHYLIP");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < testPHYLIP.size(); j++) {

            if (strcmp(SupMorfi::getIdentifier(testPHYLIP.at(j).getDetails()).c_str(),
                       compareSetEBI.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(testPHYLIP.at(j).getGene().c_str(), compareSetEBI.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 4);
    EXPECT_EQ(geneMatches, 4);
    EXPECT_EQ(testPHYLIP.size(), 4);

    std::cout << "Function Test - test_file_2_obj - PHYLIP : Complete" << std::endl;
}

//Test type PIR data extraction
TEST(function_test, test_file_2_obj_PIR) {
    //Data for the GCG, EMBL, GenBank, MSF, PHYLIP, PIR and Clustral data formats taken from the EBI website
    std::vector<GeneData> compareSetEBI = getEBIDataTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> testPIR = morfiConvert::getInData(pathPIR, "PIR");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < testPIR.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(testPIR.at(j).getDetails()).c_str(),
                       compareSetEBI.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(testPIR.at(j).getGene().c_str(), compareSetEBI.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 4);
    EXPECT_EQ(geneMatches, 4);
    EXPECT_EQ(testPIR.size(), 4);

    std::cout << "Function Test - test_file_2_obj - PIR : Complete" << std::endl;
}

//Test type Clustral data extraction
TEST(function_test, test_file_2_obj_Clustral) {
    //Data for the GCG, EMBL, GenBank, MSF, PHYLIP, PIR and Clustral data formats taken from the EBI website
    std::vector<GeneData> compareSetEBI = getEBIDataTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> testClus = morfiConvert::getInData(pathClus, "Clustral");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < testClus.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(testClus.at(j).getDetails()).c_str(),
                       compareSetEBI.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(testClus.at(j).getGene().c_str(), compareSetEBI.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 4);
    EXPECT_EQ(geneMatches, 4);
    EXPECT_EQ(testClus.size(), 4);

    std::cout << "Function Test - test_file_2_obj - Clustral : Complete" << std::endl;
}

//Test type FASTA-Report data extraction
TEST(function_test, test_file_2_obj_FR) {
    //Data for the FASTA Report data formats taken from the EBI website
    std::vector<GeneData> compareSetEBIFR = getEBIDataFASTARepTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> testFR = morfiConvert::getInData(pathFR, "FASTA-Report");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < testFR.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(testFR.at(j).getDetails()).c_str(),
                       compareSetEBIFR.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(testFR.at(j).getGene().c_str(), compareSetEBIFR.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 2);
    EXPECT_EQ(geneMatches, 2);
    EXPECT_EQ(testFR.size(), 50);

    std::cout << "Function Test - test_file_2_obj - FASTA Report : Complete" << std::endl;
}

//Test type Image MSA data extraction
TEST(function_test, test_img_2_obj) {
    int geneMatches = 0;

    std::vector<std::string> headers;
    headers.push_back("1");
    headers.push_back("2");
    headers.push_back("3");
    headers.push_back("4");
    headers.push_back("5");
    headers.push_back("6");
    headers.push_back("7");
    headers.push_back("8");
    headers.push_back("9");
    headers.push_back("10");
    headers.push_back("11");
    headers.push_back("12");
    headers.push_back("13");
    headers.push_back("14");
    headers.push_back("15");

    std::vector<std::string> geneSeqs;
    geneSeqs.push_back(
            "MA------SVSATMISTSFMPRKPAVTSL-KPIPNVGE--ALFGLKS-A--NGGKVTCMASYKVKLITPDGPIEFDCPDNVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneSeqs.push_back(
            "MA------SISGTMISTSFLPRKPAVTSL-KAISNVGE--ALFGLKS-G--RNGRITCMASYKVKLITPEGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGSVDQS");
    geneSeqs.push_back(
            "MA------SISGTMISTSFLPRKPVVTSL-KAISNVGE--ALFGLKS-G--RNGRITCMASYKVKLITPDGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGTVDQS");
    geneSeqs.push_back(
            "MAAT--TAALSGATMSTAFAPK--TPPMTAALPTNVGR--ALFGLKS-SASR-GRVTAMAAYKVTLVTPEGKQELECPDDVYILDAAEEAGIDLPYSCRAGSCSSCAGKVTSGSVNQD");
    geneSeqs.push_back(
            "MAAT--TTTMMG--MATTFVPKPQAPPMMAALPSNTGR--SLFGLKT-GSR--GGRMTMAAYKVTLVTPTGNVEFQCPDDVYILDAAEEEGIDLPYSCRAGSCSSCAGKLKTGSLNQD");
    geneSeqs.push_back(
            "MAST----ALSSAIVGTSFIRRSPAPISLRSLPSANTQ--SLFGLKS-GTARGGRVTAMATYKVKFITPEGELEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAST----ALSSAIVSTSFLRRQQTPISLRSLPFANTQ--SLFGLKS-STARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneSeqs.push_back(
            "MATT---PALYGTAVSTSFLRTQPMPMSV-TTTKAFSN--GFLGLKT-SLKRGDLAVAMASYKVKLVTPDGTQEFECPSDVYILDHAEEVGIDLPYSCRAGSCSSCAGKVVGGEVDQS");
    geneSeqs.push_back(
            "MATT---PALYGTAVSTSFMRRQPVPMSV-ATTTTTKAFPSGFGLKSVSTKRGDLAVAMATYKVKLITPEGPQEFDCPDDVYILDHAEEVGIELPYSCRAGSCSSCAGKVVNGNVNQE");
    geneSeqs.push_back(
            "MATVLGSPRAPAFFFSSSSLRAAPAPTAV--ALPAAKV--GIMGRSA-SSRR--RLRAQATYNVKLITPEGEVELQVPDDVYILDQAEEDGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAAT---------ALSMSILR---APPPCFSSPLRLRV--AVAKPLA-APMRRQLLRAQATYNVKLITPEGEVELQVPDDVYILDFAEEEGIDLPFSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAST----ALSSAIVSTSFLRRQQTPISLRSLPFANTQ--SLFGLKS-STARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneSeqs.push_back(
            "-----------------------------------------------------------ASYKVKLITPDGPIEFDCPDDVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneSeqs.push_back(
            "-----------------------------------------------------------ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "-----------------------------------------------------------ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGFVDQS");

    std::vector<std::string> geneComp;
    geneComp.push_back(
            "MASVSATMISTSFMPRKPAVTSLKPIPNVGEALFGLKSANGGKVTCMASYKVKLITPDGPIEFDCPDNVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneComp.push_back(
            "MASISGTMISTSFLPRKPAVTSLKAISNVGEALFGLKSGRNGRITCMASYKVKLITPEGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGSVDQS");
    geneComp.push_back(
            "MASISGTMISTSFLPRKPVVTSLKAISNVGEALFGLKSGRNGRITCMASYKVKLITPDGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGTVDQS");
    geneComp.push_back(
            "MAATTAALSGATMSTAFAPKTPPMTAALPTNVGRALFGLKSSASRGRVTAMAAYKVTLVTPEGKQELECPDDVYILDAAEEAGIDLPYSCRAGSCSSCAGKVTSGSVNQD");
    geneComp.push_back(
            "MAATTTTMMGMATTFVPKPQAPPMMAALPSNTGRSLFGLKTGSRGGRMTMAAYKVTLVTPTGNVEFQCPDDVYILDAAEEEGIDLPYSCRAGSCSSCAGKLKTGSLNQD");
    geneComp.push_back(
            "MASTALSSAIVGTSFIRRSPAPISLRSLPSANTQSLFGLKSGTARGGRVTAMATYKVKFITPEGELEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back(
            "MASTALSSAIVSTSFLRRQQTPISLRSLPFANTQSLFGLKSSTARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneComp.push_back(
            "MATTPALYGTAVSTSFLRTQPMPMSVTTTKAFSNGFLGLKTSLKRGDLAVAMASYKVKLVTPDGTQEFECPSDVYILDHAEEVGIDLPYSCRAGSCSSCAGKVVGGEVDQS");
    geneComp.push_back(
            "MATTPALYGTAVSTSFMRRQPVPMSVATTTTTKAFPSGFGLKSVSTKRGDLAVAMATYKVKLITPEGPQEFDCPDDVYILDHAEEVGIELPYSCRAGSCSSCAGKVVNGNVNQE");
    geneComp.push_back(
            "MATVLGSPRAPAFFFSSSSLRAAPAPTAVALPAAKVGIMGRSASSRRRLRAQATYNVKLITPEGEVELQVPDDVYILDQAEEDGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back(
            "MAATALSMSILRAPPPCFSSPLRLRVAVAKPLAAPMRRQLLRAQATYNVKLITPEGEVELQVPDDVYILDFAEEEGIDLPFSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back(
            "MASTALSSAIVSTSFLRRQQTPISLRSLPFANTQSLFGLKSSTARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneComp.push_back("ASYKVKLITPDGPIEFDCPDDVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneComp.push_back("ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back("ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGFVDQS");

    morfiConvert::removeImgTrainData();
    morfiConvert::trainImageData(geneSeqs, msaImage);

    std::vector<GeneData> test = morfiConvert::getInImageData(headers, msaImage);

    for (int i = 0; i < geneComp.size(); i++) {
        if (strcmp(geneComp.at(i).c_str(), test.at(i).getGene().c_str()) == 0) {
            geneMatches += 1;
        } else {
            std::cout << i << std::endl;
            std::cout << geneComp.at(i).c_str() << std::endl;
            std::cout << test.at(i).getGene().c_str() << std::endl;
        }
    }

    EXPECT_EQ(geneMatches, 15);

    morfiConvert::removeImgTrainData();

    std::cout << "Function Test - test_file_2_obj - MSA Image : Complete" << std::endl;

}

std::vector<GeneData> getGeneDataNCBIBLASTTextInvTest() {
    std::vector<GeneData> outRes;
    GeneData r1("AAD44166.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVILGLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGXIENY");
    outRes.push_back(r1);
    GeneData r3("YP_626379.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSNSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGMIENY");
    outRes.push_back(r3);
    GeneData r5("BAA25017.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEHPYIIIGQMASILYFSIILAFLPIAGMIENY");
    outRes.push_back(r5);

    return outRes;
}

//Test NCBI BLAST data extraction in .txt file with some invalid data
TEST(function_test, test_file_2_obj_NCBI_txt_invalid) {
    //NCBI BLAST data with come data being invalid in the original .txt file
    std::vector<GeneData> compareSetInvText = getGeneDataNCBIBLASTTextInvTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> test10 = morfiConvert::getInData(pathNCBITXTInvalid, "NCBI-BLAST");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < test10.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(test10.at(j).getDetails()).c_str(),
                       compareSetInvText.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(test10.at(j).getGene().c_str(), compareSetInvText.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 3);
    EXPECT_EQ(geneMatches, 3);
    EXPECT_EQ(test10.size(), 98);

    std::cout << "Function Test - test_file_2_obj - BLAST .txt some invalid : Complete" << std::endl;

}

std::vector<GeneData> getGeneDataNCBIBLASTOtherInvTest() {
    std::vector<GeneData> outRes;
//    GeneData r2 ("BAA25010.1","LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGMIENY");
//    outRes.push_back(r2);
    GeneData r3("YP_626379.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSNSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGMIENY");
    outRes.push_back(r3);
    GeneData r4("AAD44169.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEHPYIIIGQMASILYFSIILAFLPIAGMIENY");
    outRes.push_back(r4);
    GeneData r5("BAA25017.1",
                "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNYMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSILILGLMPFLHTSKHRSMMLRPLSQVLFWTLTMDLLTLTWIGSQPVEHPYIIIGQMASILYFSIILAFLPIAGMIENY");
    outRes.push_back(r5);

    return outRes;
}

//Test NCBI BLAST data extraction in .xml file with some invalid data
TEST(function_test, test_file_2_obj_NCBI_xml_invalid) {
    //NCBI BLAST data with come data being invalid in the original .xml and .json file
    std::vector<GeneData> compareSetInvOther = getGeneDataNCBIBLASTOtherInvTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> test11 = morfiConvert::getInData(pathNCBIXMLInvalid, "NCBI-BLAST");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < test11.size(); j++) {
            if (strcmp(SupMorfi::getIdentifier(test11.at(j).getDetails()).c_str(),
                       compareSetInvOther.at(i).getDetails().c_str()) == 0) {
                desMatches += 1;
                if (strcmp(test11.at(j).getGene().c_str(), compareSetInvOther.at(i).getGene().c_str()) == 0) {
                    geneMatches += 1;
                }
            }
        }
    }
    EXPECT_EQ(desMatches, 3);
    EXPECT_EQ(geneMatches, 3);
    EXPECT_EQ(test11.size(), 98);

    std::cout << "Function Test - test_file_2_obj - BLAST .xml some invalid : Complete" << std::endl;

}

//Test NCBI BLAST data extraction in .json file with some invalid data
TEST(function_test, test_file_2_obj_NCBI_json_invalid) {
    //NCBI BLAST data with come data being invalid in the original .xml and .json file
    std::vector<GeneData> compareSetInvOther = getGeneDataNCBIBLASTOtherInvTest();

    int geneMatches = 0;
    int desMatches = 0;

    std::vector<GeneData> test12 = morfiConvert::getInData(pathNCBIJSONInvalid, "NCBI-BLAST");
    for (int i = 1; i < 3; i++) {
        EXPECT_EQ(test12.at(i).getGene(), compareSetInvOther.at(i).getGene());
        EXPECT_EQ(SupMorfi::getIdentifier(test12.at(i).getDetails()), compareSetInvOther.at(i).getDetails());
    }
    EXPECT_EQ(test12.size(), 98);

    std::cout << "Function Test - test_file_2_obj - BLAST .json some invalid : Complete" << std::endl;

}

//Wrong file given NCBI BALST
TEST(function_test, test_file_2_obj_wrong_file_format_NCBI_BLAST) {
    try {
        std::vector<GeneData> test13 = morfiConvert::getInData(pathCusMSA, "NCBI-BLAST");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    }
    catch (std::exception &e) {
        if (strcmp(e.what(), "Incorrect input file data format") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Incorrect data format NCBI BLAST : Complete" << std::endl;
}

//Wrong file given FASTA
TEST(function_test, test_file_2_obj_wrong_file_format_FASTA) {
    try {
        std::vector<GeneData> test13 = morfiConvert::getInData(pathCusMSA, "FASTA");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    }
    catch (std::exception &e) {
        if (strcmp(e.what(), "Incorrect input file data format") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Incorrect data format FASTA : Complete" << std::endl;
}

//Wrong file given GCG
TEST(function_test, test_file_2_obj_wrong_file_format_GCG) {
    try {
        std::vector<GeneData> test13 = morfiConvert::getInData(pathCusMSA, "GCG");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    }
    catch (std::exception &e) {
        if (strcmp(e.what(), "Incorrect input file data format") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Incorrect data format GCG : Complete" << std::endl;
}

//Wrong file given EMBL
TEST(function_test, test_file_2_obj_wrong_file_format_EMBL) {
    try {
        std::vector<GeneData> test13 = morfiConvert::getInData(pathCusMSA, "EMBL");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    }
    catch (std::exception &e) {
        if (strcmp(e.what(), "Incorrect input file data format") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Incorrect data format EMBL : Complete" << std::endl;
}

//Wrong file given GenBank
TEST(function_test, test_file_2_obj_wrong_file_format_GenBank) {
    try {
        std::vector<GeneData> test13 = morfiConvert::getInData(pathCusMSA, "GenBank");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    }
    catch (std::exception &e) {
        if (strcmp(e.what(), "Incorrect input file data format") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Incorrect data format GenBank : Complete" << std::endl;
}

//Wrong file given PIR
TEST(function_test, test_file_2_obj_wrong_file_format_PIR) {
    try {
        std::vector<GeneData> test13 = morfiConvert::getInData(pathCusMSA, "PIR");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    }
    catch (std::exception &e) {
        if (strcmp(e.what(), "Incorrect input file data format") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Incorrect data format PIR : Complete" << std::endl;
}

//Wrong file given FASTA Report
TEST(function_test, test_file_2_obj_wrong_file_format_FR) {
    try {
        std::vector<GeneData> test13 = morfiConvert::getInData(pathCusMSA, "FASTA-Report");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    }
    catch (std::exception &e) {
        if (strcmp(e.what(), "Incorrect input file data format") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Incorrect data format FASTA-Report : Complete" << std::endl;
}

//Wrong path given NCBI BLAST
TEST(function_test, test_file_2_obj_wrong_path_NCBI_BLAST) {
    try {
        std::vector<GeneData> test14 = morfiConvert::getInData(wrongpath, "NCBI-BLAST");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Wrong path NCBI BLAST : Complete" << std::endl;
}

//Wrong path given MSA
TEST(function_test, test_file_2_obj_wrong_path_MSA) {
    try {
        std::vector<GeneData> test15 = morfiConvert::getInData(wrongpath, "MSA");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Wrong path MSA : Complete" << std::endl;
}

//Wrong path given FASTA
TEST(function_test, test_file_2_obj_wrong_path_FASTA) {
    try {
        std::vector<GeneData> test16 = morfiConvert::getInData(wrongpath, "FASTA");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_file_2_obj - Wrong path FASTA : Complete" << std::endl;
}

//Wrong path given BLAST
TEST(function_test, test_file_2_obj_wrong_path_BLAST) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "BLAST");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path Custom BLAST : Complete" << std::endl;
}

//Wrong path given GCG
TEST(function_test, test_file_2_obj_wrong_path_GCG) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "GCG");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path GCG : Complete" << std::endl;
}

//Wrong path given EMBL
TEST(function_test, test_file_2_obj_wrong_path_EMBL) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "EMBL");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path EMBL : Complete" << std::endl;
}

//Wrong path given PIR
TEST(function_test, test_file_2_obj_wrong_path_PIR) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "PIR");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path PIR : Complete" << std::endl;
}

//Wrong path given GenBank
TEST(function_test, test_file_2_obj_wrong_path_GenBank) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "GenBank");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path GenBank : Complete" << std::endl;
}

//Wrong path given PHYLIP
TEST(function_test, test_file_2_obj_wrong_path_PHYLIP) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "PHYLIP");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path PHYLIP : Complete" << std::endl;
}

//Wrong path given Clustral
TEST(function_test, test_file_2_obj_wrong_path_Clustral) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "Clustral");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path Clustral : Complete" << std::endl;
}

//Wrong path given MSF
TEST(function_test, test_file_2_obj_wrong_path_MSF) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "MSF");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path MSF : Complete" << std::endl;
}

//Wrong path given FASTA-Report
TEST(function_test, test_file_2_obj_wrong_path_FASTA_Report) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "FASTA-Report");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path FASTA-Report : Complete" << std::endl;
}

//Wrong path given FIND
TEST(function_test, test_file_2_obj_wrong_path_FIND) {
    try {
        std::vector<GeneData> test17 = morfiConvert::getInData(wrongpath, "FIND");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_file_2_obj - Wrong path FIND : Complete" << std::endl;
}

//ID and Description separation
TEST(function_test, test_ID_n_des_separate) {
    std::vector<std::string> test1 = SupMorfi::separateIdentifier("BAA20512.1");
    EXPECT_EQ(test1.at(0), "BAA20512.1");
    EXPECT_EQ(test1.at(1), "");

    std::vector<std::string> test2 = SupMorfi::separateIdentifier(
            "BAA20512.1 Cyprinus carpio (common carp) alpha-globin ");
    EXPECT_EQ(test2.at(0), "BAA20512.1");
    EXPECT_EQ(test2.at(1), "Cyprinus carpio (common carp) alpha-globin");

    std::cout << "Function Test - test_ID_n_des_separate : Complete" << std::endl;
}

std::vector<std::string> getOutGeneDataFASTA() {
    std::vector<std::string> out;

    std::string r1 = ">AAD44166.1 cytochrome b, partial (mitochondrion) [Elephas maximus maximus]";
    r1 += "\n";
    r1 += "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVD";
    r1 += "\n";
    r1 += "KATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLG";
    r1 += "\n";
    r1 += "DPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVILGLMPFLHTSKHRSMMLRPLSQALFWTLTMD";
    r1 += "\n";
    r1 += "LLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGXIENY";

    out.push_back(r1);

    return out;
}

//Data format conversion to FASTA
TEST(function_test, test_obj_2_obj_FASTA) {
    std::vector<std::string> comGenes = getOutGeneDataFASTA();

    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    std::vector<std::string> objOuts1 = morfiConvert::giveOutData(test1, "FASTA");
    bool matched = false;
    for (std::string compGene:objOuts1) {
        if (compGene == comGenes[0]) {
            matched = true;
        }
    }
    EXPECT_EQ(matched, true);
    EXPECT_EQ(objOuts1.size(), 100);

    std::cout << "Function Test - test_obj_2_obj - FASTA : Complete" << std::endl;
}

std::vector<std::string> getOutGeneDataIDs() {
    std::vector<std::string> out;

    out.push_back("AAD44166.1");
    out.push_back("BAA25010.1");
    out.push_back("YP_626379.1");
    out.push_back("AAD44169.1");
    out.push_back("BAA25017.1");
    out.push_back("AAD44163.1");
    out.push_back("O47885.1");
    out.push_back("NP_009291.1");
    out.push_back("AFX95036.1");
    out.push_back("AAR06164.1");

    return out;
}

//Data format conversion to NCBI IDs
TEST(function_test, test_obj_2_obj_NCBI_IDs) {
    std::vector<std::string> comIDs = getOutGeneDataIDs();

    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    std::vector<std::string> objOuts2 = morfiConvert::giveOutData(test1, "NCBI-IDs");
    int matchedIDs = 0;
    for (std::string compTranIDs:objOuts2) {
        for (std::string compMatchID:comIDs) {
            if (compMatchID == compTranIDs) {
                matchedIDs += 1;
            }
        }
    }
    EXPECT_EQ(matchedIDs, 10);
    EXPECT_EQ(objOuts2.size(), 100);

    std::cout << "Function Test - test_obj_2_obj - NCBI IDs : Complete" << std::endl;
}

std::vector<std::string> getOutGeneDataGCG() {
    std::vector<std::string> out;

    std::string r1 = "!!AA_SEQUENCE 1.0";
    r1 += "\n\n";
    r1 += "cytochrome b, partial (mitochondrion) [Elephas maximus maximus]";
    r1 += "\n\n";
    r1 += "AAD44166.1  Length: 284  Type: P  Check: 7089 ..";
    r1 += "\n\n";
    r1 += "   1 LCLYTHIGRN IYYGSYLYSE TWNTGIMLLL ITMATAFMGY VLPWGQMSFW";
    r1 += "\n\n";
    r1 += "  51 GATVITNLFS AIPYIGTNLV EWIWGGFSVD KATLNRFFAF HFILPFTMVA";
    r1 += "\n\n";
    r1 += " 101 LAGVHLTFLH ETGSNNPLGL TSDSDKIPFH PYYTIKDFLG LLILILLLLL";
    r1 += "\n\n";
    r1 += " 151 LALLSPDMLG DPDNHMPADP LNTPLHIKPE WYFLFAYAIL RSVPNKLGGV";
    r1 += "\n\n";
    r1 += " 201 LALFLSIVIL GLMPFLHTSK HRSMMLRPLS QALFWTLTMD LLTLTWIGSQ";
    r1 += "\n\n";
    r1 += " 251 PVEYPYTIIG QMASILYFSI ILAFLPIAGX IENY";
    r1 += "\n";

    out.push_back(r1);

    return out;
}

//Data format conversion to GCG
TEST(function_test, test_obj_2_obj_GCG) {
    std::vector<std::string> comGenes = getOutGeneDataGCG();

    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    std::vector<std::string> objOuts1 = morfiConvert::giveOutData(test1, "GCG");
    bool matched = false;
    for (std::string compGene:objOuts1) {
        if (compGene == comGenes[0]) {
            matched = true;
        }
    }
    EXPECT_EQ(matched, true);
    EXPECT_EQ(objOuts1.size(), 100);

    std::cout << "Function Test - test_obj_2_obj - GCG : Complete" << std::endl;
}

std::vector<std::string> getOutGeneDataEMBL() {
    std::vector<std::string> out;

    std::string r1 = "ID BAA20512.1; SV 1; linear; unassigned DNA; STD; UNC; 432 BP.";
    r1 += "\nXX\n";
    r1 += "DE Cyprinus carpio (common carp) alpha-globin";
    r1 += "\nXX\nXX\n";
    r1 += "SQ Sequence 432 BP; 97 A; 122 C; 111 G; 102 T; 0 other;";
    r1 += "\n";
    r1 += " ATGAGTCTCT CTGATAAGGA CAAGGCTGCT GTGAAAGCCC TATGGGCTAA GATCAGCCCC 60";
    r1 += "\n";
    r1 += " AAAGCCGATG ATATCGGCGC TGAAGCTCTC GGCAGAATGC TGACCGTCTA CCCTCAGACC 120";
    r1 += "\n";
    r1 += " AAGACCTACT TCGCTCACTG GGATGACCTG AGCCCTGGGT CCGGTCCTGT GAAGAAGCAT 180";
    r1 += "\n";
    r1 += " GGCAAGGTTA TCATGGGTGC AGTGGCCGAT GCCGTTTCAA AAATAGACGA CCTTGTGGGA 240";
    r1 += "\n";
    r1 += " GGTCTGGCCT CCCTGAGCGA ACTTCATGCT TCCAAGCTGC GTGTTGACCC GGCCAACTTC 300";
    r1 += "\n";
    r1 += " AAGATCCTCG CACACAATGT CATCGTGGTC ATCGGCATGC TCTTCCCTGG AGACTTCCCC 360";
    r1 += "\n";
    r1 += " CCAGAGGTTC ACATGTCAGT TGACAAGTTT TTCCAGAACT TGGCTCTGGC TCTCTCTGAG 420";
    r1 += "\n";
    r1 += " AAGTACCGCT AA 432";
    r1 += "\n//";

    out.push_back(r1);

    return out;
}

//Data format conversion to EMBL
TEST(function_test, test_obj_2_obj_EMBL) {
    std::vector<std::string> comGenes = getOutGeneDataEMBL();

    std::vector<GeneData> test1 = morfiConvert::getInData(pathGCG, "GCG");

    std::vector<std::string> objOuts1 = morfiConvert::giveOutData(test1, "EMBL");
    bool matched = false;
    for (std::string compGene:objOuts1) {
        if (compGene == comGenes[0]) {
            matched = true;
        }
    }
    EXPECT_EQ(matched, true);
    EXPECT_EQ(objOuts1.size(), 4);

    std::cout << "Function Test - test_obj_2_obj - EMBL : Complete" << std::endl;
}

std::vector<std::string> getOutGeneDataGenBank() {
    std::vector<std::string> out;

    std::string r1 = "LOCUS AAD44166.1 284 bp DNA linear UNC";
    r1 += "\n";
    r1 += "DEFINITION cytochrome b, partial (mitochondrion) [Elephas maximus maximus] .";
    r1 += "\n";
    r1 += "ACCESSION AAD44166";
    r1 += "\n";
    r1 += "ORIGIN";
    r1 += "\n";
    r1 += " 1 LCLYTHIGRN IYYGSYLYSE TWNTGIMLLL ITMATAFMGY VLPWGQMSFW GATVITNLFS";
    r1 += "\n";
    r1 += " 61 AIPYIGTNLV EWIWGGFSVD KATLNRFFAF HFILPFTMVA LAGVHLTFLH ETGSNNPLGL";
    r1 += "\n";
    r1 += " 121 TSDSDKIPFH PYYTIKDFLG LLILILLLLL LALLSPDMLG DPDNHMPADP LNTPLHIKPE";
    r1 += "\n";
    r1 += " 181 WYFLFAYAIL RSVPNKLGGV LALFLSIVIL GLMPFLHTSK HRSMMLRPLS QALFWTLTMD";
    r1 += "\n";
    r1 += " 241 LLTLTWIGSQ PVEYPYTIIG QMASILYFSI ILAFLPIAGX IENY";
    r1 += "\n//";

    out.push_back(r1);

    return out;
}

//Data format conversion to GenBank
TEST(function_test, test_obj_2_obj_Genbank) {
    std::vector<std::string> comGenes = getOutGeneDataGenBank();

    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    std::vector<std::string> objOuts1 = morfiConvert::giveOutData(test1, "GenBank");

    bool matched = false;
    for (std::string compGene:objOuts1) {
        if (compGene == comGenes[0]) {
            matched = true;
        }
    }
    EXPECT_EQ(matched, true);
    EXPECT_EQ(objOuts1.size(), 100);

    std::cout << "Function Test - test_obj_2_obj - GenBank : Complete" << std::endl;
}

std::vector<std::string> getOutGeneDataPIR() {
    std::vector<std::string> out;

    std::string r1 = ">P1;AAD44166.1";
    r1 += "\n";
    r1 += "cytochrome b, partial (mitochondrion) [Elephas maximus maximus] , 284 bases";
    r1 += "\n";
    r1 += " LCLYTHIGRN IYYGSYLYSE TWNTGIMLLL ITMATAFMGY VLPWGQMSFW";
    r1 += "\n";
    r1 += " GATVITNLFS AIPYIGTNLV EWIWGGFSVD KATLNRFFAF HFILPFTMVA";
    r1 += "\n";
    r1 += " LAGVHLTFLH ETGSNNPLGL TSDSDKIPFH PYYTIKDFLG LLILILLLLL";
    r1 += "\n";
    r1 += " LALLSPDMLG DPDNHMPADP LNTPLHIKPE WYFLFAYAIL RSVPNKLGGV";
    r1 += "\n";
    r1 += " LALFLSIVIL GLMPFLHTSK HRSMMLRPLS QALFWTLTMD LLTLTWIGSQ";
    r1 += "\n";
    r1 += " PVEYPYTIIG QMASILYFSI ILAFLPIAGX IENY*";
    r1 += "\n";

    out.push_back(r1);

    return out;
}

//Data format conversion to PIR
TEST(function_test, test_obj_2_obj_PIR) {
    std::vector<std::string> comGenes = getOutGeneDataPIR();

    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    std::vector<std::string> objOuts1 = morfiConvert::giveOutData(test1, "PIR");

    bool matched = false;
    for (std::string compGene:objOuts1) {
        if (compGene == comGenes[0]) {
            matched = true;
        }
    }
    EXPECT_EQ(matched, true);
    EXPECT_EQ(objOuts1.size(), 100);

    std::cout << "Function Test - test_obj_2_obj - PIR : Complete" << std::endl;
}

std::vector<std::string> getOutGeneDataPHYLIP() {
    std::vector<std::string> out;

    std::string r1 = "AAD44166.1LCLYTHIGRN IYYGSYLYSE TWNTGIMLLL ITMATAFMGY VLPWGQMSFW";
    r1 += "\n";
    r1 += "          GATVITNLFS AIPYIGTNLV EWIWGGFSVD KATLNRFFAF HFILPFTMVA";
    r1 += "\n";
    r1 += "          LAGVHLTFLH ETGSNNPLGL TSDSDKIPFH PYYTIKDFLG LLILILLLLL";
    r1 += "\n";
    r1 += "          LALLSPDMLG DPDNHMPADP LNTPLHIKPE WYFLFAYAIL RSVPNKLGGV";
    r1 += "\n";
    r1 += "          LALFLSIVIL GLMPFLHTSK HRSMMLRPLS QALFWTLTMD LLTLTWIGSQ";
    r1 += "\n";
    r1 += "          PVEYPYTIIG QMASILYFSI ILAFLPIAGX IENY";

    out.push_back(r1);

    return out;
}

//Data format conversion to PIR
TEST(function_test, test_obj_2_obj_PHYLIP) {
    std::vector<std::string> comGenes = getOutGeneDataPHYLIP();

    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    std::vector<std::string> objOuts1 = morfiConvert::giveOutData(test1, "PHYLIP");

    bool matched = false;
    for (std::string compGene:objOuts1) {
        if (compGene == comGenes[0]) {
            matched = true;
        }
    }
    EXPECT_EQ(matched, true);
    EXPECT_EQ(objOuts1.size(), 100);

    std::cout << "Function Test - test_obj_2_obj - PHYLIP : Complete" << std::endl;
}

//Protien data to EMBL conversion
TEST(function_test, test_obj_2_obj_EMBL_Protien) {
    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    try {
        std::vector<std::string> objOuts1 = morfiConvert::giveOutData(test1, "EMBL");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Invalid data to the output data format type") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_obj_2_obj - EMBL Protien data error : Complete" << std::endl;
}

//Invalid output data type
TEST(function_test, test_obj_2_obj_invalid_type) {
    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    try {
        std::vector<std::string> objOuts3 = morfiConvert::giveOutData(test1, "NCBI_ABC");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Unsupported output data format type") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_obj_2_obj - Invalid output data type : Complete" << std::endl;

}

//Test file write
TEST(function_test, test_obj_2_file_valid) {
    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    morfiConvert::writeOutData(test1, "res.txt", "FASTA");
    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open("res.txt");
    std::vector<std::string> sVec;

    //Temporary storage string
    std::string sTemp;

    // Add all the content in the file to sVec (To incorporate parallelism when processing the data)
    while (getline(infile, sTemp)) {
        sVec.push_back(sTemp);
    }
    //Close and clear input file stream once done
    infile.close();
    infile.clear();
    remove("res.txt");

    int linecnt = 0;
    std::vector<std::string> objOuts1 = morfiConvert::giveOutData(test1, "FASTA");
    std::vector<std::string> testObjList;
    for (std::string gene:objOuts1) {
        std::stringstream streamGene(gene);
        while (std::getline(streamGene, sTemp)) {
            for (int i = 0; i < sVec.size(); i++) {
                if (sVec[i][0] == '>') {
                    if (sVec[i] == sTemp) {
                        linecnt += 1;

                        for (int j = i + 1; j < sVec.size(); j++) {
                            if (sVec[j][0] == '>') {
                                break;
                            } else {
                                std::getline(streamGene, sTemp);
                                if (sVec[j] == sTemp) {
                                    linecnt += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    EXPECT_EQ(sVec.size(), linecnt);

    std::cout << "Function Test - test_obj_2_file - Valid file write : Complete" << std::endl;
}

//Test file write with invalid out type
TEST(function_test, test_obj_2_file_invalid_out_type) {
    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    try {
        morfiConvert::writeOutData(test1, "res.txt", "FASA");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Unsupported output data format type") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }

    std::cout << "Function Test - test_obj_2_file - Invalid output type : Complete" << std::endl;
}

//Test file write with invalid out type
TEST(function_test, test_obj_2_file_invalid_out_location) {
    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    try {
        morfiConvert::writeOutData(test1, "F:\\Achad\\res.txt", "FASTA");
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening output file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_obj_2_file - Invalid output file : Complete" << std::endl;
}

//Test file read and write
TEST(function_test, test_file_2_file) {
    std::vector<GeneData> test1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    morfiConvert::writeOutData(test1, "res1.txt", "FASTA");
    //Input file stream
    std::ifstream infile1;

    //Open file
    infile1.open("res1.txt");
    std::vector<std::string> sVec1;

    //Temporary storage string
    std::string sTemp;

    // Add all the content in the file to sVec (To incorporate parallelism when processing the data)
    while (getline(infile1, sTemp)) {
        sVec1.push_back(sTemp);
    }
    //Close and clear input file stream once done
    infile1.close();
    infile1.clear();
    remove("res1.txt");

    morfiConvert::fileConvert(pathNCBITXTValid, "res2.txt", "NCBI-BLAST", "FASTA");
    //Input file stream
    std::ifstream infile2;
    //Open file
    infile2.open("res2.txt");
    std::vector<std::string> sVec2;

    // Add all the content in the file to sVec (To incorporate parallelism when processing the data)
    while (getline(infile2, sTemp)) {
        sVec2.push_back(sTemp);
    }
    //Close and clear input file stream once done
    infile2.close();
    infile2.clear();
    remove("res2.txt");

    int linecnt = 0;
    for (int k = 0; k < sVec1.size(); k++) {
        for (int i = 0; i < sVec2.size(); i++) {
            if (sVec2[i][0] == '>') {
                if (sVec2[i] == sVec2[k]) {
                    linecnt += 1;
                    for (int j = i + 1; j < sVec2.size(); j++) {
                        if (sVec2[j][0] == '>') {
                            break;
                        } else {
                            k++;
                            if (sVec2[j] == sVec2[k]) {
                                linecnt += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    remove("res1.txt");
    remove("res2.txt");
    EXPECT_EQ(sVec2.size(), linecnt);

    std::cout << "Function Test - test_file_2_file - File read and write : Complete" << std::endl;
}

//Test identify
TEST(function_test, test_identify_valid) {
    EXPECT_EQ(morfiConvert::identify(pathNCBITXTValid), "NCBI-BLAST");
    EXPECT_EQ(morfiConvert::identify(pathNCBIJSONValid), "NCBI-BLAST");
    EXPECT_EQ(morfiConvert::identify(pathNCBIXMLValid), "NCBI-BLAST");
    EXPECT_EQ(morfiConvert::identify(pathGCG), "GCG");
    EXPECT_EQ(morfiConvert::identify(pathFASTA), "FASTA");
    EXPECT_EQ(morfiConvert::identify(pathEMBL), "EMBL");
    EXPECT_EQ(morfiConvert::identify(pathGenBank), "GenBank");
    EXPECT_EQ(morfiConvert::identify(pathPIR), "PIR");
    EXPECT_EQ(morfiConvert::identify(pathPHYLIP), "PHYLIP");
    EXPECT_EQ(morfiConvert::identify(pathClus), "Clustral");
    EXPECT_EQ(morfiConvert::identify(pathFR), "FASTA-Report");
    EXPECT_EQ(morfiConvert::identify(pathCusMSA), "Unknown");

    std::cout << "Function Test - test_identify - All valid types : Complete" << std::endl;

}

//Test identify invalid path
TEST(function_test, test_identify_invalid_path) {
    try {
        morfiConvert::identify(wrongpath);
        EXPECT_EQ("Expect exception thrown", "Exception not thrown.");
    } catch (std::exception &e) {
        if (strcmp(e.what(), "Error in opening input file. Please check if the given path is correct.") == 0) {
            EXPECT_EQ("Exception thrown", "Exception thrown");
        } else {
            EXPECT_EQ("Expect exception thrown", "Wrong exception thrown.");
        }
    }
    std::cout << "Function Test - test_identify - Wrong path : Complete" << std::endl;
}

//Test identify invalid file
TEST(function_test, test_identify_invalid_file) {
    EXPECT_EQ(morfiConvert::identify(wrongfile), "Unknown");

    std::cout << "Function Test - test_identify - Wrong file : Complete" << std::endl;
}

std::vector<GeneData> getGeneDataNCBIBLASTSegmentTest() {
    std::vector<GeneData> outRes;
    GeneData r1("AAD44166.1",
                "YTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHL");
    outRes.push_back(r1);
    GeneData r2("BAA25010.1",
                "YTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHL");
    outRes.push_back(r2);
    GeneData r3("YP_626379.1",
                "YTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHL");
    outRes.push_back(r3);
    GeneData r4("AAD44169.1",
                "YTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHL");
    outRes.push_back(r4);
    GeneData r5("BAA25017.1",
                "YTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHL");
    outRes.push_back(r5);

    return outRes;
}

//Index segment
TEST(function_test, test_index_segment) {
    std::vector<GeneData> tst = getGeneDataNCBIBLASTTest();

    std::vector<GeneData> comp = getGeneDataNCBIBLASTSegmentTest();

    std::vector<GeneData> res1 = morfiConvert::multiIndexFragment(tst, 3, 106);

    for (int i = 0; i < comp.size(); i++) {
        EXPECT_EQ(comp.at(i).getDetails(), res1.at(i).getDetails());
        EXPECT_EQ(comp.at(i).getGene(), res1.at(i).getGene());
    }

    GeneData res2 = morfiConvert::indexFragment(tst.at(0), 3, 106);
    EXPECT_EQ(comp.at(0).getDetails(), res2.getDetails());
    EXPECT_EQ(comp.at(0).getGene(), res2.getGene());

    std::cout << "Function Test - test_index_segment : Complete" << std::endl;
}

//String segment
TEST(function_test, test_string_segment) {
    std::vector<GeneData> tst = getGeneDataNCBIBLASTTest();

    std::vector<GeneData> comp = getGeneDataNCBIBLASTSegmentTest();

    std::vector<GeneData> res1 = morfiConvert::multiStringFragment(tst, "YTHIGRNI", "LAGVHLT");

    for (int i = 0; i < comp.size(); i++) {
        EXPECT_EQ(comp.at(i).getDetails(), res1.at(i).getDetails());
        EXPECT_EQ(comp.at(i).getGene(), res1.at(i).getGene());
    }

    GeneData res2 = morfiConvert::stringFragment(tst.at(0), "YTHIGRNI", "LAGVHL");
    EXPECT_EQ(comp.at(0).getDetails(), res2.getDetails());
    EXPECT_EQ(comp.at(0).getGene(), res2.getGene());

    GeneData res3 = morfiConvert::stringFragment(tst.at(1), "YTHIGRNI", "LAGVHL");
    EXPECT_EQ(comp.at(1).getDetails(), res3.getDetails());
    EXPECT_EQ(comp.at(1).getGene(), res3.getGene());

    std::cout << "Function Test - test_string_segment : Complete" << std::endl;
}

//**********************************************************************************************************************
//Performance Test
//**********************************************************************************************************************

//Local store
TEST(performance_test, test_local_store) {
    remove("points.morfi");
    long origSize = 0;
    //1.9 GB

    std::ifstream origfile;
    origfile.open(datapath, std::ifstream::ate | std::ifstream::binary);
    origSize = origfile.tellg();

    clock_t t1, t2;
    t1 = clock();
    morfiConvert::configLocal(datapath);
    t2 = clock();
    float diff((float) t2 - (float) t1);
    std::cout << "1.9GB conversion time :" << diff / 1000 << "s" << std::endl;
    EXPECT_TRUE(diff < 1000);

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open("points.morfi");

    std::string s;
    std::string sGene;

    int j = 1;

    std::string geneFound = "";

    long newSizes = 0;
    if (infile.good()) {
        while (getline(infile, s)) {
            //String without the final pointer value
            std::string shaveLastPoint = s.substr(0, s.find_last_of(" "));

            //Get the file path of the local database file
            std::string dbPath = shaveLastPoint.substr(0, shaveLastPoint.find_last_of(" "));

            std::ifstream datafile;
            datafile.open(dbPath);

            std::cout << dbPath << " : " << datafile.tellg() << std::endl;
            newSizes += datafile.tellg();

            datafile.close();
            remove(dbPath.c_str());
        }
        remove("points.morfi");
    } else {
        EXPECT_EQ(true, false);
    }
    if (newSizes < origSize) {
        EXPECT_EQ("New sizes are smaller", "New sizes are smaller");
    } else {
        EXPECT_EQ("Original is smaller", "New sizes are smaller");
    }

    std::cout << "Function Test - test_local_store : Complete" << std::endl;

}

//Retrieve data
TEST(performance_test, test_local_retrieve) {
    remove("points.morfi");
    //1.9 GB
    morfiConvert::configLocal(datapath);


    clock_t t1, t2;
    t1 = clock();
    std::string test1 = morfiConvert::fullGene("ARO40711.1", "all");
    t2 = clock();
    float diff((float) t2 - (float) t1);
    std::cout << "1.9GB data retrieival time :" << diff / 1000 << "s" << std::endl;

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open("points.morfi");

    std::string s;

    long newSizes = 0;
    if (infile.good()) {
        while (getline(infile, s)) {
            //String without the final pointer value
            std::string shaveLastPoint = s.substr(0, s.find_last_of(" "));

            //Get the file path of the local database file
            std::string dbPath = shaveLastPoint.substr(0, shaveLastPoint.find_last_of(" "));

            remove(dbPath.c_str());
        }
        remove("points.morfi");
    } else {
        EXPECT_EQ(true, false);
    }
    if (diff < 1000) {
        EXPECT_EQ("Retrieval time is less than 1s", "Retrieval time is less than 1s");
    } else {
        EXPECT_EQ("Retrieval time is more than 1s", "Retrieval time is less than 1s");
    }

    std::cout << "Function Test - test_local_retrieve : Complete" << std::endl;
}

//Test file 2 obj conversions
TEST(performance_test, test_file_2_obj) {
    int minThreadNum = 0;
    int minTime = 13000;

    for (int threadNum = 1; threadNum <= 10; threadNum += 1) {
        std::cout << "For " + std::to_string(threadNum) + " threads." << std::endl;
        clock_t t1, t2;
        t1 = clock();
        std::vector<GeneData> res1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST", threadNum);
        t2 = clock();
        float diff1((float) t2 - (float) t1);
        std::cout << "Time to convert a NCBI-BLAST txt file :" << diff1 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff1 < 1000);

        t1 = clock();
        std::vector<GeneData> res2 = morfiConvert::getInData(pathNCBIXMLValid, "NCBI-BLAST", threadNum);
        t2 = clock();
        float diff2((float) t2 - (float) t1);
        std::cout << "Time to convert a NCBI-BLAST xml file :" << diff2 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff2 < 1000);

        t1 = clock();
        std::vector<GeneData> res3 = morfiConvert::getInData(pathNCBIJSONValid, "NCBI-BLAST", threadNum);
        t2 = clock();
        float diff3((float) t2 - (float) t1);
        std::cout << "Time to convert a NCBI-BLAST json file :" << diff3 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff3 < 1000);

        t1 = clock();
        std::vector<GeneData> res4 = morfiConvert::getInData(pathCusMSA, "MSA");
        t2 = clock();
        float diff4((float) t2 - (float) t1);
        std::cout << "Time to convert a MSA file :" << diff4 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff4 < 1000);

        t1 = clock();
        std::vector<GeneData> res5 = morfiConvert::getInData(pathCusBLAST, "BLAST", threadNum);
        t2 = clock();
        float diff5((float) t2 - (float) t1);
        std::cout << "Time to convert a Custom BLAST file :" << diff5 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff5 < 1000);

        t1 = clock();
        std::vector<GeneData> res6 = morfiConvert::getInData(pathFASTA, "FASTA", threadNum);
        t2 = clock();
        float diff6((float) t2 - (float) t1);
        std::cout << "Time to convert a FASTA file :" << diff6 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff6 < 1000);

        t1 = clock();
        std::vector<GeneData> res7 = morfiConvert::getInData(pathGCG, "GCG", threadNum);
        t2 = clock();
        float diff7((float) t2 - (float) t1);
        std::cout << "Time to convert a GCG file :" << diff7 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff7 < 1000);

        t1 = clock();
        std::vector<GeneData> res8 = morfiConvert::getInData(pathEMBL, "EMBL", threadNum);
        t2 = clock();
        float diff8((float) t2 - (float) t1);
        std::cout << "Time to convert a EMBL file :" << diff8 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff8 < 1000);

        t1 = clock();
        std::vector<GeneData> res9 = morfiConvert::getInData(pathGenBank, "GenBank", threadNum);
        t2 = clock();
        float diff9((float) t2 - (float) t1);
        std::cout << "Time to convert a GenBank file :" << diff9 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff9 < 1000);

        t1 = clock();
        std::vector<GeneData> res10 = morfiConvert::getInData(pathPIR, "PIR", threadNum);
        t2 = clock();
        float diff10((float) t2 - (float) t1);
        std::cout << "Time to convert a PIR file :" << diff10 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff10 < 1000);

        t1 = clock();
        std::vector<GeneData> res11 = morfiConvert::getInData(pathPHYLIP, "PHYLIP", threadNum);
        t2 = clock();
        float diff11((float) t2 - (float) t1);
        std::cout << "Time to convert a PHYLIP file :" << diff11 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff11 < 1000);

        t1 = clock();
        std::vector<GeneData> res12 = morfiConvert::getInData(pathClus, "Clustral", threadNum);
        t2 = clock();
        float diff12((float) t2 - (float) t1);
        std::cout << "Time to convert a Clustral file :" << diff12 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff12 < 1000);

        t1 = clock();
        std::vector<GeneData> res13 = morfiConvert::getInData(pathMSF, "MSF", threadNum);
        t2 = clock();
        float diff13((float) t2 - (float) t1);
        std::cout << "Time to convert a MSF file :" << diff13 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff13 < 1000);

        t1 = clock();
        std::vector<GeneData> res14 = morfiConvert::getInData(pathFR, "FASTA-Report", threadNum);
        t2 = clock();
        float diff14((float) t2 - (float) t1);
        std::cout << "Time to convert a FASTA-Report file :" << diff14 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff14 < 1000);

        int sumTime = diff1 + diff2 + diff3 + diff4 + diff5 + diff6 + diff7 + diff8 + diff9 + diff10 + diff11 + diff12 +
                      diff13;
        if (minTime >= (sumTime)) {
            minTime = sumTime;
            minThreadNum = threadNum;
        }
    }
    std::cout << minThreadNum << " threads produce the best result of " << minTime << "ms" << std::endl;

    std::cout << "Function Test - test_file_2_obj : Complete" << std::endl;
}

//Test obj 2 obj conversions
TEST(performance_test, test_obj_2_obj) {
    std::vector<GeneData> res1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");
    std::vector<GeneData> res2 = morfiConvert::getInData(pathGCG, "GCG");

    int minThreadNum = 0;
    int minTime = 6040;

    for (int threadNum = 1; threadNum <= 10; threadNum += 1) {
        std::cout << "For " + std::to_string(threadNum) + " threads." << std::endl;
        clock_t t1, t2;
        t1 = clock();
        std::vector<std::string> hld1 = morfiConvert::giveOutData(res1, "FASTA", threadNum);
        t2 = clock();
        float diff1((float) t2 - (float) t1);
        std::cout << "Time to convert GeneData objs to FASTA objs :" << diff1 << "ms" << std::endl;
        EXPECT_TRUE(diff1 < 1000);

        t1 = clock();
        std::vector<std::string> hld2 = morfiConvert::giveOutData(res1, "NCBI-IDs", threadNum);
        t2 = clock();
        float diff2((float) t2 - (float) t1);
        std::cout << "Time to convert GeneData objs to NCBI-IDs :" << diff2 << "ms" << std::endl;
        EXPECT_TRUE(diff2 < 1000);

        t1 = clock();
        std::vector<std::string> hld3 = morfiConvert::giveOutData(res1, "GCG", threadNum);
        t2 = clock();
        float diff3((float) t2 - (float) t1);
        std::cout << "Time to convert GeneData objs to GCG :" << diff3 << "ms" << std::endl;
        EXPECT_TRUE(diff3 < 1000);

        t1 = clock();
        std::vector<std::string> hld4 = morfiConvert::giveOutData(res2, "EMBL", threadNum);
        t2 = clock();
        float diff4((float) t2 - (float) t1);
        std::cout << "Time to convert GeneData objs to EMBL :" << diff4 << "ms" << std::endl;
        EXPECT_TRUE(diff4 < 40);

        t1 = clock();
        std::vector<std::string> hld5 = morfiConvert::giveOutData(res1, "GenBank", threadNum);
        t2 = clock();
        float diff5((float) t2 - (float) t1);
        std::cout << "Time to convert GeneData objs to GenBank :" << diff5 << "ms" << std::endl;
        EXPECT_TRUE(diff5 < 1000);

        t1 = clock();
        std::vector<std::string> hld6 = morfiConvert::giveOutData(res1, "PIR", threadNum);
        t2 = clock();
        float diff6((float) t2 - (float) t1);
        std::cout << "Time to convert GeneData objs to PIR :" << diff6 << "ms" << std::endl;
        EXPECT_TRUE(diff6 < 1000);

        t1 = clock();
        std::vector<std::string> hld7 = morfiConvert::giveOutData(res1, "PHYLIP", threadNum);
        t2 = clock();
        float diff7((float) t2 - (float) t1);
        std::cout << "Time to convert GeneData objs to PHYLIP :" << diff7 << "ms" << std::endl;
        EXPECT_TRUE(diff7 < 1000);

        int sumTime = diff1 + diff2 + diff3 + diff4 + diff5 + diff6 + diff7;
        if (minTime >= (sumTime)) {
            minTime = sumTime;
            minThreadNum = threadNum;
        }
    }
    std::cout << minThreadNum << " threads produce the best result of " << minTime << "ms" << std::endl;

    std::cout << "Function Test - test_obj_2_obj : Complete" << std::endl;
}

//Test obj 2 file conversions
TEST(performance_test, test_obj_2_file) {
    std::vector<GeneData> res1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");

    clock_t t1, t2;
    t1 = clock();
    morfiConvert::writeOutData(res1, "res1.txt", "FASTA");
    t2 = clock();
    remove("res1.txt");
    float diff1((float) t2 - (float) t1);
    std::cout << "Time to convert GeneData objs to FASTA file :" << diff1 << "ms" << std::endl;
    EXPECT_TRUE(diff1 < 1000);

    std::cout << "Function Test - test_obj_2_file : Complete" << std::endl;

}

//Test file type identification
TEST(performance_test, test_identify) {
    clock_t t1, t2;
    t1 = clock();
    std::string res1 = morfiConvert::identify(pathNCBITXTValid);
    t2 = clock();
    float diff1((float) t2 - (float) t1);
    std::cout << "Time to Identify FASTA txt file :" << diff1 << "ms" << std::endl;
    EXPECT_TRUE(diff1 < 1000);

    t1 = clock();
    std::string res2 = morfiConvert::identify(pathNCBIXMLValid);
    t2 = clock();
    float diff2((float) t2 - (float) t1);
    std::cout << "Time to Identify FASTA xml file :" << diff2 << "ms" << std::endl;
    EXPECT_TRUE(diff2 < 1000);

    t1 = clock();
    std::string res3 = morfiConvert::identify(pathNCBIJSONValid);
    t2 = clock();
    float diff3((float) t2 - (float) t1);
    std::cout << "Time to Identify FASTA json file :" << diff3 << "ms" << std::endl;
    EXPECT_TRUE(diff3 < 1000);

    t1 = clock();
    std::string res4 = morfiConvert::identify(pathGCG);
    t2 = clock();
    float diff4((float) t2 - (float) t1);
    std::cout << "Time to Identify GCG file :" << diff4 << "ms" << std::endl;
    EXPECT_TRUE(diff4 < 1000);

    t1 = clock();
    std::string res5 = morfiConvert::identify(pathFASTA);
    t2 = clock();
    float diff5((float) t2 - (float) t1);
    std::cout << "Time to Identify FASTA file :" << diff5 << "ms" << std::endl;
    EXPECT_TRUE(diff5 < 1000);

    t1 = clock();
    std::string res6 = morfiConvert::identify(pathEMBL);
    t2 = clock();
    float diff6((float) t2 - (float) t1);
    std::cout << "Time to Identify EMBL file :" << diff6 << "ms" << std::endl;
    EXPECT_TRUE(diff6 < 1000);

    t1 = clock();
    std::string res7 = morfiConvert::identify(pathGenBank);
    t2 = clock();
    float diff7((float) t2 - (float) t1);
    std::cout << "Time to Identify GenBank file :" << diff7 << "ms" << std::endl;
    EXPECT_TRUE(diff7 < 1000);

    t1 = clock();
    std::string res8 = morfiConvert::identify(pathPIR);
    t2 = clock();
    float diff8((float) t2 - (float) t1);
    std::cout << "Time to Identify PIR file :" << diff8 << "ms" << std::endl;
    EXPECT_TRUE(diff8 < 1000);

    t1 = clock();
    std::string res9 = morfiConvert::identify(pathPHYLIP);
    t2 = clock();
    float diff9((float) t2 - (float) t1);
    std::cout << "Time to Identify PHYLIP file :" << diff9 << "ms" << std::endl;
    EXPECT_TRUE(diff9 < 1000);

    t1 = clock();
    std::string res10 = morfiConvert::identify(pathClus);
    t2 = clock();
    float diff10((float) t2 - (float) t1);
    std::cout << "Time to Identify Clustral file :" << diff10 << "ms" << std::endl;
    EXPECT_TRUE(diff10 < 1000);

    t1 = clock();
    std::string res11 = morfiConvert::identify(pathFR);
    t2 = clock();
    float diff11((float) t2 - (float) t1);
    std::cout << "Time to Identify FASTA-Report file :" << diff11 << "ms" << std::endl;
    EXPECT_TRUE(diff11 < 1000);

    std::cout << "Function Test - test_identify : Complete" << std::endl;
}

//Test image character training
TEST(performance_test, test_img_2_obj_train) {
    clock_t t1, t2;

    std::vector<std::string> headers;
    headers.push_back("1");
    headers.push_back("2");
    headers.push_back("3");
    headers.push_back("4");
    headers.push_back("5");
    headers.push_back("6");
    headers.push_back("7");
    headers.push_back("8");
    headers.push_back("9");
    headers.push_back("10");
    headers.push_back("11");
    headers.push_back("12");
    headers.push_back("13");
    headers.push_back("14");
    headers.push_back("15");

    std::vector<std::string> geneSeqs;
    morfiConvert::removeImgTrainData();

    geneSeqs.push_back(
            "MA------SVSATMISTSFMPRKPAVTSL-KPIPNVGE--ALFGLKS-A--NGGKVTCMASYKVKLITPDGPIEFDCPDNVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneSeqs.push_back(
            "MA------SISGTMISTSFLPRKPAVTSL-KAISNVGE--ALFGLKS-G--RNGRITCMASYKVKLITPEGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGSVDQS");
    geneSeqs.push_back(
            "MA------SISGTMISTSFLPRKPVVTSL-KAISNVGE--ALFGLKS-G--RNGRITCMASYKVKLITPDGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGTVDQS");
    geneSeqs.push_back(
            "MAAT--TAALSGATMSTAFAPK--TPPMTAALPTNVGR--ALFGLKS-SASR-GRVTAMAAYKVTLVTPEGKQELECPDDVYILDAAEEAGIDLPYSCRAGSCSSCAGKVTSGSVNQD");
    geneSeqs.push_back(
            "MAAT--TTTMMG--MATTFVPKPQAPPMMAALPSNTGR--SLFGLKT-GSR--GGRMTMAAYKVTLVTPTGNVEFQCPDDVYILDAAEEEGIDLPYSCRAGSCSSCAGKLKTGSLNQD");
    geneSeqs.push_back(
            "MAST----ALSSAIVGTSFIRRSPAPISLRSLPSANTQ--SLFGLKS-GTARGGRVTAMATYKVKFITPEGELEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAST----ALSSAIVSTSFLRRQQTPISLRSLPFANTQ--SLFGLKS-STARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneSeqs.push_back(
            "MATT---PALYGTAVSTSFLRTQPMPMSV-TTTKAFSN--GFLGLKT-SLKRGDLAVAMASYKVKLVTPDGTQEFECPSDVYILDHAEEVGIDLPYSCRAGSCSSCAGKVVGGEVDQS");
    geneSeqs.push_back(
            "MATT---PALYGTAVSTSFMRRQPVPMSV-ATTTTTKAFPSGFGLKSVSTKRGDLAVAMATYKVKLITPEGPQEFDCPDDVYILDHAEEVGIELPYSCRAGSCSSCAGKVVNGNVNQE");
    geneSeqs.push_back(
            "MATVLGSPRAPAFFFSSSSLRAAPAPTAV--ALPAAKV--GIMGRSA-SSRR--RLRAQATYNVKLITPEGEVELQVPDDVYILDQAEEDGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAAT---------ALSMSILR---APPPCFSSPLRLRV--AVAKPLA-APMRRQLLRAQATYNVKLITPEGEVELQVPDDVYILDFAEEEGIDLPFSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAST----ALSSAIVSTSFLRRQQTPISLRSLPFANTQ--SLFGLKS-STARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneSeqs.push_back(
            "-----------------------------------------------------------ASYKVKLITPDGPIEFDCPDDVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneSeqs.push_back(
            "-----------------------------------------------------------ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "-----------------------------------------------------------ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGFVDQS");


    t1 = clock();
    morfiConvert::trainImageData(geneSeqs, msaImage);
    t2 = clock();
    float diff1((float) t2 - (float) t1);
    std::cout << "Time to train characters in image :" << diff1 << "ms" << std::endl;
    EXPECT_TRUE(diff1 < 6000);

    std::cout << "Function Test - test_img_2_obj_train : Complete" << std::endl;
}

//Test image to file conversion
TEST(performance_test, test_img_2_obj) {
    clock_t t1, t2;

    std::vector<std::string> headers;
    headers.push_back("1");
    headers.push_back("2");
    headers.push_back("3");
    headers.push_back("4");
    headers.push_back("5");
    headers.push_back("6");
    headers.push_back("7");
    headers.push_back("8");
    headers.push_back("9");
    headers.push_back("10");
    headers.push_back("11");
    headers.push_back("12");
    headers.push_back("13");
    headers.push_back("14");
    headers.push_back("15");

    std::vector<std::string> geneSeqs;
    geneSeqs.push_back(
            "MA------SVSATMISTSFMPRKPAVTSL-KPIPNVGE--ALFGLKS-A--NGGKVTCMASYKVKLITPDGPIEFDCPDNVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneSeqs.push_back(
            "MA------SISGTMISTSFLPRKPAVTSL-KAISNVGE--ALFGLKS-G--RNGRITCMASYKVKLITPEGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGSVDQS");
    geneSeqs.push_back(
            "MA------SISGTMISTSFLPRKPVVTSL-KAISNVGE--ALFGLKS-G--RNGRITCMASYKVKLITPDGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGTVDQS");
    geneSeqs.push_back(
            "MAAT--TAALSGATMSTAFAPK--TPPMTAALPTNVGR--ALFGLKS-SASR-GRVTAMAAYKVTLVTPEGKQELECPDDVYILDAAEEAGIDLPYSCRAGSCSSCAGKVTSGSVNQD");
    geneSeqs.push_back(
            "MAAT--TTTMMG--MATTFVPKPQAPPMMAALPSNTGR--SLFGLKT-GSR--GGRMTMAAYKVTLVTPTGNVEFQCPDDVYILDAAEEEGIDLPYSCRAGSCSSCAGKLKTGSLNQD");
    geneSeqs.push_back(
            "MAST----ALSSAIVGTSFIRRSPAPISLRSLPSANTQ--SLFGLKS-GTARGGRVTAMATYKVKFITPEGELEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAST----ALSSAIVSTSFLRRQQTPISLRSLPFANTQ--SLFGLKS-STARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneSeqs.push_back(
            "MATT---PALYGTAVSTSFLRTQPMPMSV-TTTKAFSN--GFLGLKT-SLKRGDLAVAMASYKVKLVTPDGTQEFECPSDVYILDHAEEVGIDLPYSCRAGSCSSCAGKVVGGEVDQS");
    geneSeqs.push_back(
            "MATT---PALYGTAVSTSFMRRQPVPMSV-ATTTTTKAFPSGFGLKSVSTKRGDLAVAMATYKVKLITPEGPQEFDCPDDVYILDHAEEVGIELPYSCRAGSCSSCAGKVVNGNVNQE");
    geneSeqs.push_back(
            "MATVLGSPRAPAFFFSSSSLRAAPAPTAV--ALPAAKV--GIMGRSA-SSRR--RLRAQATYNVKLITPEGEVELQVPDDVYILDQAEEDGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAAT---------ALSMSILR---APPPCFSSPLRLRV--AVAKPLA-APMRRQLLRAQATYNVKLITPEGEVELQVPDDVYILDFAEEEGIDLPFSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAST----ALSSAIVSTSFLRRQQTPISLRSLPFANTQ--SLFGLKS-STARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneSeqs.push_back(
            "-----------------------------------------------------------ASYKVKLITPDGPIEFDCPDDVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneSeqs.push_back(
            "-----------------------------------------------------------ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "-----------------------------------------------------------ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGFVDQS");

    std::vector<std::string> geneComp;
    geneComp.push_back(
            "MASVSATMISTSFMPRKPAVTSLKPIPNVGEALFGLKSANGGKVTCMASYKVKLITPDGPIEFDCPDNVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneComp.push_back(
            "MASISGTMISTSFLPRKPAVTSLKAISNVGEALFGLKSGRNGRITCMASYKVKLITPEGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGSVDQS");
    geneComp.push_back(
            "MASISGTMISTSFLPRKPVVTSLKAISNVGEALFGLKSGRNGRITCMASYKVKLITPDGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGTVDQS");
    geneComp.push_back(
            "MAATTAALSGATMSTAFAPKTPPMTAALPTNVGRALFGLKSSASRGRVTAMAAYKVTLVTPEGKQELECPDDVYILDAAEEAGIDLPYSCRAGSCSSCAGKVTSGSVNQD");
    geneComp.push_back(
            "MAATTTTMMGMATTFVPKPQAPPMMAALPSNTGRSLFGLKTGSRGGRMTMAAYKVTLVTPTGNVEFQCPDDVYILDAAEEEGIDLPYSCRAGSCSSCAGKLKTGSLNQD");
    geneComp.push_back(
            "MASTALSSAIVGTSFIRRSPAPISLRSLPSANTQSLFGLKSGTARGGRVTAMATYKVKFITPEGELEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back(
            "MASTALSSAIVSTSFLRRQQTPISLRSLPFANTQSLFGLKSSTARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneComp.push_back(
            "MATTPALYGTAVSTSFLRTQPMPMSVTTTKAFSNGFLGLKTSLKRGDLAVAMASYKVKLVTPDGTQEFECPSDVYILDHAEEVGIDLPYSCRAGSCSSCAGKVVGGEVDQS");
    geneComp.push_back(
            "MATTPALYGTAVSTSFMRRQPVPMSVATTTTTKAFPSGFGLKSVSTKRGDLAVAMATYKVKLITPEGPQEFDCPDDVYILDHAEEVGIELPYSCRAGSCSSCAGKVVNGNVNQE");
    geneComp.push_back(
            "MATVLGSPRAPAFFFSSSSLRAAPAPTAVALPAAKVGIMGRSASSRRRLRAQATYNVKLITPEGEVELQVPDDVYILDQAEEDGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back(
            "MAATALSMSILRAPPPCFSSPLRLRVAVAKPLAAPMRRQLLRAQATYNVKLITPEGEVELQVPDDVYILDFAEEEGIDLPFSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back(
            "MASTALSSAIVSTSFLRRQQTPISLRSLPFANTQSLFGLKSSTARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneComp.push_back("ASYKVKLITPDGPIEFDCPDDVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneComp.push_back("ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back("ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGFVDQS");

    morfiConvert::removeImgTrainData();
    morfiConvert::trainImageData(geneSeqs, msaImage);


    t1 = clock();
    std::vector<GeneData> test = morfiConvert::getInImageData(headers, msaImage);
    t2 = clock();
    morfiConvert::removeImgTrainData();
    float diff1((float) t2 - (float) t1);
    std::cout << "Time to Identify data in image :" << diff1 << "ms" << std::endl;
    EXPECT_TRUE(diff1 < 4000);

    std::cout << "Function Test - test_img_2_obj : Complete" << std::endl;
}

//**********************************************************************************************************************
//Load test
//**********************************************************************************************************************

//Retrieve data
TEST(load_test, test_local_retrieve) {
    remove("points.morfi");
    //1.9 GB
    morfiConvert::configLocal(datapath);


    clock_t t1, t2;
    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string test1 = morfiConvert::fullGene("ARO40711.1", "all");
    }
    t2 = clock();
    float diff((float) t2 - (float) t1);
    EXPECT_TRUE(diff < 60000);
    std::cout << "1.9GB 100x data retrieival time :" << diff / 1000 << "s" << std::endl;

    //Input file stream
    std::ifstream infile;

    //Open file
    infile.open("points.morfi");

    std::string s;

    long newSizes = 0;
    if (infile.good()) {
        while (getline(infile, s)) {
            //String without the final pointer value
            std::string shaveLastPoint = s.substr(0, s.find_last_of(" "));

            //Get the file path of the local database file
            std::string dbPath = shaveLastPoint.substr(0, shaveLastPoint.find_last_of(" "));

            remove(dbPath.c_str());
        }
        remove("points.morfi");
    } else {
        EXPECT_EQ(true, false);
    }
    if (diff < 1000) {
        EXPECT_EQ("Retrieval time is less than 1s", "Retrieval time is less than 1s");
    } else {
        EXPECT_EQ("Retrieval time is more than 1s", "Retrieval time is less than 1s");
    }

    std::cout << "Load Test - test_local_retrieve : Complete" << std::endl;
}

//File 2 obj conversion
TEST(load_test, test_file_2_obj) {
    int minThreadNum = 0;
    int minTime = 60000;

    for (int threadNum = 1; threadNum <= 10; threadNum += 1) {
        std::cout << "For " + std::to_string(threadNum) + " threads." << std::endl;
        clock_t t1, t2;
        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST", threadNum);
        }
        t2 = clock();
        float diff1((float) t2 - (float) t1);
        std::cout << "Time to convert 100 NCBI-BLAST txt files :" << diff1 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff1 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res2 = morfiConvert::getInData(pathNCBIXMLValid, "NCBI-BLAST", threadNum);
        }
        t2 = clock();
        float diff2((float) t2 - (float) t1);
        std::cout << "Time to convert 100 NCBI-BLAST xml files :" << diff2 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff2 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res3 = morfiConvert::getInData(pathNCBIJSONValid, "NCBI-BLAST", threadNum);
        }
        t2 = clock();
        float diff3((float) t2 - (float) t1);
        std::cout << "Time to convert 100 NCBI-BLAST json files :" << diff3 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff3 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res4 = morfiConvert::getInData(pathCusMSA, "MSA");
        }
        t2 = clock();
        float diff4((float) t2 - (float) t1);
        std::cout << "Time to convert 100 MSA files :" << diff4 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff4 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res5 = morfiConvert::getInData(pathCusBLAST, "BLAST", threadNum);
        }
        t2 = clock();
        float diff5((float) t2 - (float) t1);
        std::cout << "Time to convert 100 Custom BLAST files :" << diff5 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff5 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res6 = morfiConvert::getInData(pathFASTA, "FASTA", threadNum);
        }
        t2 = clock();
        float diff6((float) t2 - (float) t1);
        std::cout << "Time to convert 100 FASTA files :" << diff6 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff6 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res7 = morfiConvert::getInData(pathGCG, "GCG", threadNum);
        }
        t2 = clock();
        float diff7((float) t2 - (float) t1);
        std::cout << "Time to convert 100 GCG files :" << diff7 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff7 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res8 = morfiConvert::getInData(pathEMBL, "EMBL", threadNum);
        }
        t2 = clock();
        float diff8((float) t2 - (float) t1);
        std::cout << "Time to convert 100 EMBL files :" << diff8 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff8 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res9 = morfiConvert::getInData(pathGenBank, "GenBank", threadNum);
        }
        t2 = clock();
        float diff9((float) t2 - (float) t1);
        std::cout << "Time to convert 100 GenBank files :" << diff9 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff9 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res10 = morfiConvert::getInData(pathPIR, "PIR", threadNum);
        }
        t2 = clock();
        float diff10((float) t2 - (float) t1);
        std::cout << "Time to convert 100 PIR files :" << diff10 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff10 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res11 = morfiConvert::getInData(pathPHYLIP, "PHYLIP", threadNum);
        }
        t2 = clock();
        float diff11((float) t2 - (float) t1);
        std::cout << "Time to convert 100 PHYLIP files :" << diff11 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff11 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res12 = morfiConvert::getInData(pathClus, "Clustral", threadNum);
        }
        t2 = clock();
        float diff12((float) t2 - (float) t1);
        std::cout << "Time to convert 100 Clustral files :" << diff12 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff12 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res13 = morfiConvert::getInData(pathMSF, "MSF", threadNum);
        }
        t2 = clock();
        float diff13((float) t2 - (float) t1);
        std::cout << "Time to convert 100 MSF files :" << diff13 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff13 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<GeneData> res14 = morfiConvert::getInData(pathFR, "FASTA-Report", threadNum);
        }
        t2 = clock();
        float diff14((float) t2 - (float) t1);
        std::cout << "Time to convert 100 FASTA-Report files :" << diff14 / 1000 << "s" << std::endl;
        EXPECT_TRUE(diff14 < 60000);

        int sumTime = diff1 + diff2 + diff3 + diff4 + diff5 + diff6 + diff7 + diff8 + diff9 + diff10 + diff11 + diff12 +
                      diff13;
        if (minTime >= (sumTime)) {
            minTime = sumTime;
            minThreadNum = threadNum;
        }
    }
    std::cout << minThreadNum << " threads produce the best result of " << minTime << "ms" << std::endl;

    std::cout << "Load Test - test_file_2_obj : Complete" << std::endl;
}

//Obj 2 obj conversions
TEST(load_test, test_obj_2_obj) {
    std::vector<GeneData> res1 = morfiConvert::getInData(pathNCBITXTValid, "NCBI-BLAST");
    std::vector<GeneData> res2 = morfiConvert::getInData(pathGCG, "GCG");

    int minThreadNum = 0;
    int minTime = 60400;

    for (int threadNum = 1; threadNum <= 10; threadNum += 1) {
        std::cout << "For " + std::to_string(threadNum) + " threads." << std::endl;

        clock_t t1, t2;
        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<std::string> hld1 = morfiConvert::giveOutData(res1, "FASTA", threadNum);
        }
        t2 = clock();
        float diff1((float) t2 - (float) t1);
        std::cout << "Time to convert 100 GeneData obj sets to FASTA objs :" << diff1 << "ms" << std::endl;
        EXPECT_TRUE(diff1 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<std::string> hld2 = morfiConvert::giveOutData(res1, "NCBI-IDs", threadNum);
        }
        t2 = clock();
        float diff2((float) t2 - (float) t1);
        std::cout << "Time to convert 100 GeneData obj sets to NCBI-IDs :" << diff2 << "ms" << std::endl;
        EXPECT_TRUE(diff2 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<std::string> hld3 = morfiConvert::giveOutData(res1, "GCG", threadNum);
        }
        t2 = clock();
        float diff3((float) t2 - (float) t1);
        std::cout << "Time to convert 100 GeneData objs to GCG :" << diff3 << "ms" << std::endl;
        EXPECT_TRUE(diff3 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<std::string> hld4 = morfiConvert::giveOutData(res2, "EMBL", threadNum);
        }
        t2 = clock();
        float diff4((float) t2 - (float) t1);
        std::cout << "Time to convert 100 GeneData objs to EMBL :" << diff4 << "ms" << std::endl;
        EXPECT_TRUE(diff4 < 2400);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<std::string> hld5 = morfiConvert::giveOutData(res1, "GenBank", threadNum);
        }
        t2 = clock();
        float diff5((float) t2 - (float) t1);
        std::cout << "Time to convert 100 GeneData objs to GenBank :" << diff5 << "ms" << std::endl;
        EXPECT_TRUE(diff5 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<std::string> hld6 = morfiConvert::giveOutData(res1, "PIR", threadNum);
        }
        t2 = clock();
        float diff6((float) t2 - (float) t1);
        std::cout << "Time to convert 100 GeneData objs to PIR :" << diff6 << "ms" << std::endl;
        EXPECT_TRUE(diff6 < 60000);

        t1 = clock();
#pragma omp parallel for
        for (int i = 0; i < 100; i++) {
            std::vector<std::string> hld7 = morfiConvert::giveOutData(res1, "PHYLIP", threadNum);
        }
        t2 = clock();
        float diff7((float) t2 - (float) t1);
        std::cout << "Time to convert 100 GeneData objs to PHYLIP :" << diff7 << "ms" << std::endl;
        EXPECT_TRUE(diff7 < 60000);

        int sumTime = diff1 + diff2 + diff3 + diff4 + diff5 + diff6 + diff7;
        if (minTime >= (sumTime)) {
            minTime = sumTime;
            minThreadNum = threadNum;
        }
    }
    std::cout << minThreadNum << " threads produce the best result of " << minTime << "ms" << std::endl;

    std::cout << "Load Test - test_obj_2_obj : Complete" << std::endl;
}

//Obj 2 file conversions
TEST(load_test, test_obj_2_file) {
    std::string path1 = "F:\\Achademic\\CS\\Semester 05\\Software Engineering Project\\Tests\\Protien\\Valid Total (NCBI)\\Test1.txt";
    std::vector<GeneData> res1 = morfiConvert::getInData(path1, "NCBI-BLAST");

    clock_t t1, t2;
    t1 = clock();
    for (int i = 0; i < 100; i++) {
        morfiConvert::writeOutData(res1, "res1.txt", "FASTA");
    }
    remove("res1.txt");
    t2 = clock();
    float diff1((float) t2 - (float) t1);
    std::cout << "Time to convert 100 GeneData obj sets to FASTA files :" << diff1 << "ms" << std::endl;
    EXPECT_TRUE(diff1 < 60000);
}

//Identify
TEST(load_test, test_identify) {
    clock_t t1, t2;
    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res1 = morfiConvert::identify(pathNCBITXTValid);
    }
    t2 = clock();
    float diff1((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 FASTA txt files :" << diff1 << "ms" << std::endl;
    EXPECT_TRUE(diff1 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res2 = morfiConvert::identify(pathNCBIXMLValid);
    }
    t2 = clock();
    float diff2((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 FASTA xml files :" << diff2 << "ms" << std::endl;
    EXPECT_TRUE(diff2 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res3 = morfiConvert::identify(pathNCBIJSONValid);
    }
    t2 = clock();
    float diff3((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 FASTA json files :" << diff3 << "ms" << std::endl;
    EXPECT_TRUE(diff3 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res4 = morfiConvert::identify(pathGCG);
    }
    t2 = clock();
    float diff4((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 GCG files :" << diff4 << "ms" << std::endl;
    EXPECT_TRUE(diff4 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res5 = morfiConvert::identify(pathFASTA);
    }
    t2 = clock();
    float diff5((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 FASTA files :" << diff5 << "ms" << std::endl;
    EXPECT_TRUE(diff5 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res6 = morfiConvert::identify(pathEMBL);
    }
    t2 = clock();
    float diff6((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 EMBL files :" << diff6 << "ms" << std::endl;
    EXPECT_TRUE(diff6 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res7 = morfiConvert::identify(pathGenBank);
    }
    t2 = clock();
    float diff7((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 GenBank files :" << diff7 << "ms" << std::endl;
    EXPECT_TRUE(diff7 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res8 = morfiConvert::identify(pathPIR);
    }
    t2 = clock();
    float diff8((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 PIR files :" << diff8 << "ms" << std::endl;
    EXPECT_TRUE(diff8 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res9 = morfiConvert::identify(pathPHYLIP);
    }
    t2 = clock();
    float diff9((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 PHYLIP files :" << diff9 << "ms" << std::endl;
    EXPECT_TRUE(diff9 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res10 = morfiConvert::identify(pathClus);
    }
    t2 = clock();
    float diff10((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 Clustral files :" << diff10 << "ms" << std::endl;
    EXPECT_TRUE(diff10 < 60000);

    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::string res11 = morfiConvert::identify(pathFR);
    }
    t2 = clock();
    float diff11((float) t2 - (float) t1);
    std::cout << "Time to Identify 100 FASTA-Report files :" << diff11 << "ms" << std::endl;
    EXPECT_TRUE(diff11 < 60000);

    std::cout << "Load Test - test_identify : Complete" << std::endl;
}

//Test image to file conversion
TEST(load_test, test_img_2_obj) {
    clock_t t1, t2;

    std::vector<std::string> headers;
    headers.push_back("1");
    headers.push_back("2");
    headers.push_back("3");
    headers.push_back("4");
    headers.push_back("5");
    headers.push_back("6");
    headers.push_back("7");
    headers.push_back("8");
    headers.push_back("9");
    headers.push_back("10");
    headers.push_back("11");
    headers.push_back("12");
    headers.push_back("13");
    headers.push_back("14");
    headers.push_back("15");

    std::vector<std::string> geneSeqs;
    geneSeqs.push_back(
            "MA------SVSATMISTSFMPRKPAVTSL-KPIPNVGE--ALFGLKS-A--NGGKVTCMASYKVKLITPDGPIEFDCPDNVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneSeqs.push_back(
            "MA------SISGTMISTSFLPRKPAVTSL-KAISNVGE--ALFGLKS-G--RNGRITCMASYKVKLITPEGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGSVDQS");
    geneSeqs.push_back(
            "MA------SISGTMISTSFLPRKPVVTSL-KAISNVGE--ALFGLKS-G--RNGRITCMASYKVKLITPDGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGTVDQS");
    geneSeqs.push_back(
            "MAAT--TAALSGATMSTAFAPK--TPPMTAALPTNVGR--ALFGLKS-SASR-GRVTAMAAYKVTLVTPEGKQELECPDDVYILDAAEEAGIDLPYSCRAGSCSSCAGKVTSGSVNQD");
    geneSeqs.push_back(
            "MAAT--TTTMMG--MATTFVPKPQAPPMMAALPSNTGR--SLFGLKT-GSR--GGRMTMAAYKVTLVTPTGNVEFQCPDDVYILDAAEEEGIDLPYSCRAGSCSSCAGKLKTGSLNQD");
    geneSeqs.push_back(
            "MAST----ALSSAIVGTSFIRRSPAPISLRSLPSANTQ--SLFGLKS-GTARGGRVTAMATYKVKFITPEGELEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAST----ALSSAIVSTSFLRRQQTPISLRSLPFANTQ--SLFGLKS-STARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneSeqs.push_back(
            "MATT---PALYGTAVSTSFLRTQPMPMSV-TTTKAFSN--GFLGLKT-SLKRGDLAVAMASYKVKLVTPDGTQEFECPSDVYILDHAEEVGIDLPYSCRAGSCSSCAGKVVGGEVDQS");
    geneSeqs.push_back(
            "MATT---PALYGTAVSTSFMRRQPVPMSV-ATTTTTKAFPSGFGLKSVSTKRGDLAVAMATYKVKLITPEGPQEFDCPDDVYILDHAEEVGIELPYSCRAGSCSSCAGKVVNGNVNQE");
    geneSeqs.push_back(
            "MATVLGSPRAPAFFFSSSSLRAAPAPTAV--ALPAAKV--GIMGRSA-SSRR--RLRAQATYNVKLITPEGEVELQVPDDVYILDQAEEDGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAAT---------ALSMSILR---APPPCFSSPLRLRV--AVAKPLA-APMRRQLLRAQATYNVKLITPEGEVELQVPDDVYILDFAEEEGIDLPFSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "MAST----ALSSAIVSTSFLRRQQTPISLRSLPFANTQ--SLFGLKS-STARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneSeqs.push_back(
            "-----------------------------------------------------------ASYKVKLITPDGPIEFDCPDDVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneSeqs.push_back(
            "-----------------------------------------------------------ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneSeqs.push_back(
            "-----------------------------------------------------------ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGFVDQS");

    std::vector<std::string> geneComp;
    geneComp.push_back(
            "MASVSATMISTSFMPRKPAVTSLKPIPNVGEALFGLKSANGGKVTCMASYKVKLITPDGPIEFDCPDNVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneComp.push_back(
            "MASISGTMISTSFLPRKPAVTSLKAISNVGEALFGLKSGRNGRITCMASYKVKLITPEGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGSVDQS");
    geneComp.push_back(
            "MASISGTMISTSFLPRKPVVTSLKAISNVGEALFGLKSGRNGRITCMASYKVKLITPDGPIEFECPDDVYILDQAEEEGHDLPYSCRAGSCSSCAGKVTAGTVDQS");
    geneComp.push_back(
            "MAATTAALSGATMSTAFAPKTPPMTAALPTNVGRALFGLKSSASRGRVTAMAAYKVTLVTPEGKQELECPDDVYILDAAEEAGIDLPYSCRAGSCSSCAGKVTSGSVNQD");
    geneComp.push_back(
            "MAATTTTMMGMATTFVPKPQAPPMMAALPSNTGRSLFGLKTGSRGGRMTMAAYKVTLVTPTGNVEFQCPDDVYILDAAEEEGIDLPYSCRAGSCSSCAGKLKTGSLNQD");
    geneComp.push_back(
            "MASTALSSAIVGTSFIRRSPAPISLRSLPSANTQSLFGLKSGTARGGRVTAMATYKVKFITPEGELEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back(
            "MASTALSSAIVSTSFLRRQQTPISLRSLPFANTQSLFGLKSSTARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneComp.push_back(
            "MATTPALYGTAVSTSFLRTQPMPMSVTTTKAFSNGFLGLKTSLKRGDLAVAMASYKVKLVTPDGTQEFECPSDVYILDHAEEVGIDLPYSCRAGSCSSCAGKVVGGEVDQS");
    geneComp.push_back(
            "MATTPALYGTAVSTSFMRRQPVPMSVATTTTTKAFPSGFGLKSVSTKRGDLAVAMATYKVKLITPEGPQEFDCPDDVYILDHAEEVGIELPYSCRAGSCSSCAGKVVNGNVNQE");
    geneComp.push_back(
            "MATVLGSPRAPAFFFSSSSLRAAPAPTAVALPAAKVGIMGRSASSRRRLRAQATYNVKLITPEGEVELQVPDDVYILDQAEEDGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back(
            "MAATALSMSILRAPPPCFSSPLRLRVAVAKPLAAPMRRQLLRAQATYNVKLITPEGEVELQVPDDVYILDFAEEEGIDLPFSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back(
            "MASTALSSAIVSTSFLRRQQTPISLRSLPFANTQSLFGLKSSTARGGRVTAMATYKVKFITPEGEQEVECEEDVYVLDAAEEAGLDLPYSCRAGSCSSCAGKVVSGSIDQS");
    geneComp.push_back("ASYKVKLITPDGPIEFDCPDDVYILDQAEEAGHDLPYSCRAGSCSSCAGKIAGGAVDQT");
    geneComp.push_back("ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGSVDQS");
    geneComp.push_back("ATYKVKFITPEGEQEVECDDDVYVLDAAEEAGIDLPYSCRAGSCSSCAGKVVSGFVDQS");

    morfiConvert::removeImgTrainData();
    morfiConvert::trainImageData(geneSeqs, msaImage);


    t1 = clock();
#pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        std::vector<GeneData> test = morfiConvert::getInImageData(headers, msaImage);
    }
    t2 = clock();
    float diff1((float) t2 - (float) t1);
    std::cout << "Time to Identify data in 100 images :" << diff1 << "ms" << std::endl;
    EXPECT_TRUE(diff1 < 240000);

    std::cout << "Load Test - test_img_2_obj : Complete" << std::endl;
}

//**********************************************************************************************************************
//Regression test
//**********************************************************************************************************************

//OpenMP
TEST(regression_test, test_openmp) {
    try {
        long sm = 0;
#pragma omp parallel for num_threads(10)
        for (int i = 0; i < 10000000; i++) {
#pragma omp critical
            {
                sm += i;
            }
        }

        EXPECT_EQ("Exception not thrown", "Exception not thrown");
    } catch (std::exception &e) {
        EXPECT_EQ("Exception not thrown", "Exception thrown");
    }

    std::cout << "Regression Test - test_openmp : Complete" << std::endl;
}

//**********************************************************************************************************************
//Initiate tests
//**********************************************************************************************************************

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
