/*
 * Copyright (c) 2025 Yao LEI.
 * All rights reserved.
 *
 * This source code is licensed under the MIT License found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>
#include <map>
#include <future>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <htslib/bgzf.h>

using namespace std;

const int BLOCK_SIZE = 100000; // Number of lines per block (tune as needed)

struct BlockResult {
    size_t block_index;
    std::vector<std::string> lines;
};

vector<string> split(const string& str, char delimiter) {
    vector<string> tokens;
    stringstream ss(str);
    string token;

    while (getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}

string replace_string(const string& str, const string& from, const string& to) {
    if (from.empty())
        return str;
    string result = str;
    size_t start_pos = 0;
    while ((start_pos = result.find(from, start_pos)) != string::npos) {
        result.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
    return result;
}

vector<string> create_CSQ_header_vector(const char* CSQheaderRow) {
    string replaced_CSQheaderRow1 = replace_string(CSQheaderRow, "\n", ""); // remove the \n in the end.
    string replaced_CSQheaderRow2 = replace_string(replaced_CSQheaderRow1, "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: ", ""); // remove the unwanted string in the row. This part need to be generalized in the future.
    string replaced_CSQheaderRow3 = replace_string(replaced_CSQheaderRow2, "\">", ""); // remove the unwanted substring in the row.
    vector<string> CSQ_header_vector = split(replaced_CSQheaderRow3, '|');

    return CSQ_header_vector;
}

void process_block_to_file(size_t block_index, const vector<string>& block_lines, const vector<string>& CSQ_header_vec, const string& temp_filename) {
    ofstream temp_out(temp_filename);
    for (const auto& buffer : block_lines) {
        vector<string> token_vec = split(buffer, '\t'); // Split the line into columns by tab

        // Check if CSQ header is available and the line has at least 8 columns
        if (CSQ_header_vec.size() != 0 && token_vec.size() >= 8) {
            string CHROM = token_vec[0]; // Extract chromosome
            string POS = token_vec[1];   // Extract position
            string REF = token_vec[3];   // Extract reference allele
            string ALT = token_vec[4];   // Extract alternate allele
            string INFO = token_vec[7];  // Extract INFO field

            // Instead of writing to file, collect output lines
            vector<string> info_vec = split(INFO, ';'); // Split INFO field by ';' to get individual info fields

            for (const auto& i : info_vec) { // Iterate over each info field
                if (i.substr(0, 4) == "CSQ=") { // If the field is the CSQ annotation
                    vector<string> CSQ_vec = split(replace_string(i, "CSQ=", ""), ','); // Remove "CSQ=" and split by ',' to get each CSQ annotation

                    for (const auto& j : CSQ_vec) { // Iterate over each CSQ annotation
                        vector<string> CSQ_sub_vec = split(j, '|'); // Split the CSQ annotation by '|' to get individual terms

                        // If the number of terms is one less than the header, append an empty string to fix it
                        if (CSQ_sub_vec.size() == CSQ_header_vec.size() - 1) {
                            CSQ_sub_vec.push_back("");
                        }

                        // Start building the output line with CHROM, POS, REF, ALT
                        string output_line = CHROM + "\t" + POS + "\t" + REF + "\t" + ALT + "\t";

                        // Add each CSQ term to the output line, replacing empty terms with "."
                        for (size_t k = 0; k < CSQ_sub_vec.size() - 1; ++k) {
                            if (CSQ_sub_vec[k] == "") CSQ_sub_vec[k] = ".";
                            output_line += CSQ_sub_vec[k] + "\t";
                        }

                        // Handle the last CSQ term, replacing empty with "."
                        if (CSQ_sub_vec[CSQ_sub_vec.size() - 1] == "") CSQ_sub_vec[CSQ_sub_vec.size() - 1] = ".";
                        output_line += CSQ_sub_vec[CSQ_sub_vec.size() - 1] + "\n";

                        // Instead of collecting output_lines, write each line directly:
                        temp_out << output_line;
                    }
                }
            }
        }
    }
    temp_out.close();
}

void split_vep_into_tsv_parallel(const char* input_vcf_gz, const char* output_prefix, const char* tmp_folder) {
    gzFile file = gzopen(input_vcf_gz, "rb"); // Open the input gzipped VCF file for reading
    if (file == NULL) { // Check if the file was opened successfully
        cerr << "Error: Failed to open file: " << input_vcf_gz << "." << endl; // Print error if not
        exit(1); // Exit the program with error
    }

    // gzFile output = gzopen((string(output_prefix) + ".tsv.gz").c_str(), "wb"); // Open the output gzipped TSV file for writing
    BGZF* output = bgzf_open((string(output_prefix) + ".tsv.gz").c_str(), "w9");  // Open the output gzipped TSV file for writing
    if (!output) { // Check if the output file was opened successfully
        cerr << "Failed to open file for writing." << endl; // Print error if not
        exit(1); // Exit the program with error
    }

    const int buffer_size = 7 * 1024 * 1024; // 1Mb buffer size for reading lines
    char buffer[buffer_size]; // Buffer to hold each line read from the file
    vector<string> CSQ_header_vec; // Vector to store the CSQ header fields
    string header_line; // String to store the output header line
    vector<future<void>> futures; // Vector to store futures for async processing
    vector<string> temp_filenames; // Vector to store temp filenames

    size_t block_index = 0; // Index to keep track of block number
    vector<string> block_lines; // Vector to store lines for the current block

    // Read the input file line by line
    while (gzgets(file, buffer, buffer_size) != NULL) {
        vector<string> token_vec = split(buffer, '\t'); // Split the line into columns by tab

        // Get the CSQ header from VCF header row
        if (token_vec[0].substr(0, 14) == "##INFO=<ID=CSQ") {
            CSQ_header_vec = create_CSQ_header_vector(buffer); // Parse the CSQ header into a vector
            header_line = "CHROM\tPOS\tREF\tALT\t"; // Start building the output header line
            for (size_t i = 0; i < CSQ_header_vec.size() - 1; ++i) {
                header_line += CSQ_header_vec[i] + "\t"; // Add each CSQ header field except the last, with tab
            }
            header_line += CSQ_header_vec[CSQ_header_vec.size() - 1] + "\n"; // Add the last CSQ header field and newline
            // gzwrite(output, header_line.c_str(), header_line.size()); // Write the header line to the output file
            ssize_t written = bgzf_write(output, header_line.c_str(), header_line.size()); // Write the header line to the output file
            continue; // Skip to the next line
        }

        // If the line is not a header line, add it to the current block
        if (token_vec[0].substr(0, 1) != "#") block_lines.push_back(buffer);

        // If the block is full, process it asynchronously
        if (block_lines.size() == BLOCK_SIZE) {
            string temp_filename = string(tmp_folder) + "/" + string(output_prefix) + ".block" + to_string(block_index) + ".tmp";
            temp_filenames.push_back(temp_filename);
            futures.push_back(async(launch::async, process_block_to_file, block_index, block_lines, CSQ_header_vec, temp_filename));
            block_lines.clear(); // Clear the block for the next set of lines
            ++block_index; // Increment the block index
        }
    }
    
    // Process any remaining lines that didn't fill a complete block
    if (!block_lines.empty()) {
        string temp_filename = string(tmp_folder) + "/" + string(output_prefix) + ".block" + to_string(block_index) + ".tmp";
        temp_filenames.push_back(temp_filename);
        futures.push_back(async(launch::async, process_block_to_file, block_index, block_lines, CSQ_header_vec, temp_filename));
    }

    // Wait for all async tasks to finish
    for (auto& f : futures) f.get();

    // Write all output lines in order to the output file
    for (const auto& temp_filename : temp_filenames) {
        ifstream temp_in(temp_filename);
        string line;
        while (getline(temp_in, line)) {
            line += "\n";
            // gzwrite(output, line.c_str(), line.size());
            ssize_t written = bgzf_write(output, line.c_str(), line.size());
        }
        temp_in.close();
        remove(temp_filename.c_str()); // Delete temp file
    }

    gzclose(file);   // Close the input file
    // gzclose(output); // Close the output file
    bgzf_close(output); // Close the output file
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <input.vcf.gz> <output_prefix> <tmp_folder>" << endl;
        return 1;
    }
    const char* input_vcf_gz = argv[1];
    const char* output_prefix = argv[2];
    const char* tmp_folder = argv[3];

    split_vep_into_tsv_parallel(input_vcf_gz, output_prefix, tmp_folder);

    return 0;
}
