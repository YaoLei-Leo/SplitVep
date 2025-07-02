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

using namespace std;

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

void write_INFO_CSQ(string INFO_string, string CHROM, string POS, string REF, string ALT, vector<string> CSQ_header_vec, gzFile output) {
    vector<string> info_vec = split(INFO_string, ';');

    for (const auto& i : info_vec) { // 1. Find the CSQ part.
        if (i.substr(0, 4) == "CSQ=") {
            vector<string> CSQ_vec = split(replace_string(i, "CSQ=", ""), ','); // Remove the unwanted substring.

            for (const auto& j : CSQ_vec) {
                vector<string> CSQ_sub_vec = split(j, '|');

                // VEP annotate use | as separator in the CSQ field. However, it the final annotation term is NA, VEP writes the CSQ end with | than space or NA. Therefore, here, if the size of INFO CSQ field is 1 term less than the CSQ header. I append a " " into the vector to fix it.
                if (CSQ_sub_vec.size() == CSQ_header_vec.size() - 1) {
                    CSQ_sub_vec.push_back("");
                }

                string output_line = CHROM + "," + POS + "," + REF + "," + ALT + ",";
                for (size_t k = 0; k < CSQ_sub_vec.size() - 1; ++k) {
                    output_line += CSQ_sub_vec[k] + ",";
                }
                output_line += CSQ_sub_vec[CSQ_sub_vec.size() - 1] + "\n";

                gzwrite(output, output_line.c_str(), output_line.size());
            }
        }
    }
}

void split_vep_into_csv(const char* input_vcf_gz, const char* output_prefix) {
    gzFile file = gzopen(input_vcf_gz, "rb");
    if (file == NULL) {
        cerr << "Error: Failed to open file: " << input_vcf_gz << "." << endl;
        exit(1);
    }

    gzFile output = gzopen((string(output_prefix) + ".csv.gz").c_str(), "wb"); // Create output gz file.
    if (!output) {
        std::cerr << "Failed to open file for writing." << endl;
        exit(1);
    }

    const int buffer_size = 1024 * 1024;
    char buffer[buffer_size]; // 1MB buffer
    vector<string> CSQ_header_vec;
    while (gzgets(file, buffer, buffer_size) != NULL) {
        vector<string> token_vec = split(buffer, '\t');
        // Create CSQ header vector.
        if (token_vec[0].substr(0, 14) == "##INFO=<ID=CSQ") { // If this row is the CSQ header.
            cout << "INFO CSQ identified." << endl;
            CSQ_header_vec = create_CSQ_header_vector(buffer); // Convert the INFO CSQ row into CSQ header vector.

            // Write the header into output file.
            string header_line = "CHROM,POS,REF,ALT,";
            for (size_t i = 0; i < CSQ_header_vec.size() - 1; ++i) {
                header_line += CSQ_header_vec[i] + ",";
            }
            header_line += CSQ_header_vec[CSQ_header_vec.size() - 1] + "\n";
            gzwrite(output, header_line.c_str(), header_line.size());
        }

        if ( token_vec[0].substr(0, 1) != "#" && CSQ_header_vec.size() != 0 && token_vec.size() >= 8) { // 8 columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.
            string CHROM = token_vec[0];
            string POS = token_vec[1];
            // string ID = token_vec[2];
            string REF = token_vec[3];
            string ALT = token_vec[4];
            // string QUAL = token_vec[5];
            // string FILTER = token_vec[6];
            string INFO = token_vec[7];

            write_INFO_CSQ(INFO, CHROM, POS, REF, ALT, CSQ_header_vec, output); // Write the INFO CSQ into output in duplicate mode.
        }
    }

    gzclose(file);
    gzclose(output);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input.vcf.gz> <output_prefix>" << endl;
        return 1;
    }
    const char* input_vcf_gz = argv[1];
    const char* output_prefix = argv[2];

    split_vep_into_csv(input_vcf_gz, output_prefix);

    return 0;
}