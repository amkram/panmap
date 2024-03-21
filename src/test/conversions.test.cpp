#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <stack>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include "PangenomeMAT.hpp"

#include <algorithm>
#include "../place.hpp"
#include "../pmi.hpp"
#include "../util.hpp"
#include "../tree.hpp"
#include "../genotype.hpp"
#include <cmath>
#include <htslib/sam.h>

#include "../conversion.hpp"

using namespace std;
using namespace tree;
namespace fs = boost::filesystem;




bool areFilesEqual(const std::string& file1_path, const std::string& file2_path) {
    std::ifstream file1(file1_path, std::ios::binary);
    std::ifstream file2(file2_path, std::ios::binary);

    if (!file1.is_open() || !file2.is_open()) {
        std::cerr << "Error: Couldn't open one or both files." << std::endl;
        return false;
    }

    char byte1, byte2;
    while (file1.get(byte1) && file2.get(byte2)) {
        if (byte1 != byte2) {
            // Files are not equal if bytes don't match
            file1.close();
            file2.close();
            cerr << "BYTES NOT EQUAL: " << (int)byte1 << " " << (int)byte2 << "\n";
            return false;
        }
    }

    file1.close();
    file2.close();
    return true;
}








BOOST_AUTO_TEST_CASE(_createBam) {

    cout << "\n\nCreate bam tests\n";

    for(int i = 0; i < 4; i++) {
        string sam_file_path = "../src/test/data/conversion_test_data/bam_input/test" + to_string(i) + "_R1.sam";


        cout << "Testing the conversion of file test" << i << "_R1.sam to a bam\n";




        //OPEN SAM
        // Open the file
        std::ifstream file(sam_file_path);
        
        // Check if the file is opened successfully
        if (!file.is_open()) {
            std::cerr << "Error: Failed to open the file." << std::endl;
            return;
        }
        
        // Vector to store lines as char pointers
        std::vector<char *> samAlignments;
        
        // String to store header lines
        std::string samHeader;
        
        // Read lines from the file
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue; // Skip empty lines
            
            // Allocate memory for the line and copy the content
            char *line_cstr = (char *)malloc(line.size() + 1); // +1 for null terminator
            std::strcpy(line_cstr, line.c_str());
            
            // Check if the line is a header line
            if (line[0] == '@') {
                // Concatenate header lines
                samHeader += line + '\n';
                free(line_cstr); // Free memory for header lines since we don't add them to the vector
            } else {
                // Add non-header lines to the vector
                samAlignments.push_back(line_cstr);
            }
        }
        
        // Close the file
        file.close();
        
        



        //PREPARE FOR createBam
        std::string bamFileName = "../src/test/data/conversion_test_data/bam_output/test" + to_string(i) + "_R1.bam";

        sam_hdr_t *header;
        bam1_t **bamRecords;




        //CALL FUNCTION
        createBam(
            samAlignments,
            samHeader,
            bamFileName,

            header,
            bamRecords
        );
        
        


        //CHECK IF WE GOT EXPECTED FILES
        std::string expected_bamFileName = "../src/test/data/conversion_test_data/bam_ex_output/test" + to_string(i) + "_R1.bam";

        BOOST_TEST(areFilesEqual(expected_bamFileName, bamFileName));

        cerr << "Finished subtest\n\n";
        
    }
    
}