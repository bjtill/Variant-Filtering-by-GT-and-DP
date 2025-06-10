//BT June 3, 2025
// g++ -o vcf_filter vcf_filter.cpp -lz
//-lz above compiles with zlib support

//Command-line options:

//-a: Allow missing data (./.) as neutral
//-m: Accept multiple alleles (e.g., 2/2, 0/2, 1/3, etc.)
//--min-depth integer min depth
//--max-depth integer maximum depth 
//-i: Input VCF file
//-s: Sample file
//-o: Output file (optional, defaults to stdout)


//# Basic filtering
//./vcf_filter -i input.vcf.gz -s samples.txt

//# Allow missing data and multiple alleles
//./vcf_filter -i input.vcf.gz -s samples.txt -a -m -o filtered.vcf

//with depth
//./vcf_filter -i test1.vcf.gz -s samples.txt -a -m --min-depth 40 --max-depth 200 -o filtered4depth.vcf

//# Show help
//./vcf_filter -h


//How it works:

//Filter Code Logic:

//When -m is used, homozygous alternative accepts any identical non-zero alleles (1/1, 2/2, 3/3)
//Heterozygous accepts any combination where one allele is 0
//Missing data (./.) is accepted as neutral only when -a flag is used


//Processing:

//Reads VCF header to map sample names to column positions
//For each variant line, checks if all specified samples match their filter codes
//Only outputs variants where ALL samples meet their criteria



//The program handles edge cases like missing samples, invalid genotype formats, and malformed VCF files with appropriate error messages.



#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <getopt.h>
#include <zlib.h>

struct Sample {
    std::string name;
    int code;
};

struct GenotypeCall {
    int allele1;
    int allele2;
    bool is_missing;
    
    GenotypeCall() : allele1(-1), allele2(-1), is_missing(true) {}
    GenotypeCall(int a1, int a2) : allele1(a1), allele2(a2), is_missing(false) {}
};

class VCFFilter {
private:
    std::vector<Sample> samples;
    bool allow_missing;
    bool multiple_alleles;
    std::map<std::string, int> sample_positions;
    int min_depth;
    int max_depth;
    
public:
    VCFFilter(const std::vector<Sample>& sample_list, bool allow_miss, bool multi_allele, 
              int min_dp = -1, int max_dp = -1)
        : samples(sample_list), allow_missing(allow_miss), multiple_alleles(multi_allele),
          min_depth(min_dp), max_depth(max_dp) {}
    
    GenotypeCall parseGenotype(const std::string& gt_str) {
        if (gt_str == "./." || gt_str == ".|.") {
            return GenotypeCall(); // missing
        }
        
        // Find separator (either / for unphased or | for phased)
        size_t sep_pos = gt_str.find('/');
        if (sep_pos == std::string::npos) {
            sep_pos = gt_str.find('|');
        }
        
        if (sep_pos == std::string::npos) {
            return GenotypeCall(); // invalid format
        }
        
        try {
            int allele1 = std::stoi(gt_str.substr(0, sep_pos));
            int allele2 = std::stoi(gt_str.substr(sep_pos + 1));
            return GenotypeCall(allele1, allele2);
        } catch (const std::exception&) {
            return GenotypeCall(); // parsing error
        }
    }
    
    bool isHomozygousReference(const GenotypeCall& gt) {
        return !gt.is_missing && gt.allele1 == 0 && gt.allele2 == 0;
    }
    
    bool isHomozygousAlternative(const GenotypeCall& gt) {
        if (gt.is_missing) return false;
        
        if (multiple_alleles) {
            return gt.allele1 == gt.allele2 && gt.allele1 > 0;
        } else {
            return gt.allele1 == 1 && gt.allele2 == 1;
        }
    }
    
    bool isHeterozygous(const GenotypeCall& gt) {
        if (gt.is_missing) return false;
        
        if (multiple_alleles) {
            // For multiple alleles, heterozygous means different alleles with at least one being 0
            return gt.allele1 != gt.allele2 && (gt.allele1 == 0 || gt.allele2 == 0);
        } else {
            // Standard case: one allele is 0, the other is 1 (accepts both 0/1|1 and 1/0|0)
            return (gt.allele1 == 0 && gt.allele2 == 1) || (gt.allele1 == 1 && gt.allele2 == 0);
        }
    }
    
    int extractDepth(const std::string& sample_data, const std::string& format_string) {
        // Parse FORMAT string to find DP position
        std::istringstream format_stream(format_string);
        std::string format_field;
        int dp_position = -1;
        int position = 0;
        
        while (std::getline(format_stream, format_field, ':')) {
            if (format_field == "DP") {
                dp_position = position;
                break;
            }
            position++;
        }
        
        if (dp_position == -1) {
            return -1; // DP not found in FORMAT
        }
        
        // Parse sample data to extract DP value
        std::istringstream sample_stream(sample_data);
        std::string field;
        position = 0;
        
        while (std::getline(sample_stream, field, ':')) {
            if (position == dp_position) {
                try {
                    if (field == "." || field.empty()) {
                        return -1; // Missing depth
                    }
                    return std::stoi(field);
                } catch (const std::exception&) {
                    return -1; // Invalid depth value
                }
            }
            position++;
        }
        
        return -1; // DP field not found in sample data
    }
    
    bool matchesCode(const GenotypeCall& gt, int code) {
        if (gt.is_missing && allow_missing) {
            return true; // missing data is neutral when -a flag is used
        }
        
        if (gt.is_missing && !allow_missing) {
            return false;
        }
        
        switch (code) {
            case 1: // homozygous reference (0/0)
                return isHomozygousReference(gt);
            case 2: // homozygous alternative (1/1, 2/2, 3/3 if -m)
                return isHomozygousAlternative(gt);
            case 3: // heterozygous (0/1, 0/2, 1/2, etc. if -m)
                return isHeterozygous(gt);
            case 4: // heterozygous or homozygous reference
                return isHeterozygous(gt) || isHomozygousReference(gt);
            case 5: // homozygous alternative or heterozygous
                return isHomozygousAlternative(gt) || isHeterozygous(gt);
            case 6: // homozygous reference or homozygous alternative
                return isHomozygousReference(gt) || isHomozygousAlternative(gt);
            case 7: // anything
                return true;
            default:
                return false;
        }
    }
    
    void setSamplePositions(const std::string& header_line) {
        std::istringstream iss(header_line);
        std::string token;
        int pos = 0;
        
        // Skip the first 9 columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
        for (int i = 0; i < 9; i++) {
            iss >> token;
        }
        
        // Map sample names to their positions
        while (iss >> token) {
            sample_positions[token] = pos;
            pos++;
        }
    }
    
    bool filterVariant(const std::string& line) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> fields;
        
        // Split the line into fields
        while (iss >> token) {
            fields.push_back(token);
        }
        
        if (fields.size() < 10) {
            return false; // Invalid VCF line
        }
        
        // Get FORMAT field (column 8, 0-indexed)
        std::string format_field = fields[8];
        
        // Extract genotype fields (starting from column 9, 0-indexed)
        std::vector<std::string> genotypes;
        for (size_t i = 9; i < fields.size(); i++) {
            genotypes.push_back(fields[i]);
        }
        
        // Collect depths for all samples in our filter list
        std::vector<int> sample_depths;
        
        // Check each sample
        for (const auto& sample : samples) {
            auto pos_it = sample_positions.find(sample.name);
            if (pos_it == sample_positions.end()) {
                std::cerr << "Warning: Sample " << sample.name << " not found in VCF header\n";
                continue;
            }
            
            int sample_pos = pos_it->second;
            if (sample_pos >= static_cast<int>(genotypes.size())) {
                std::cerr << "Warning: Sample position out of range\n";
                continue;
            }
            
            std::string sample_data = genotypes[sample_pos];
            
            // Extract GT field (assume it's the first field in FORMAT)
            size_t colon_pos = sample_data.find(':');
            std::string gt_str = (colon_pos != std::string::npos) ? 
                                sample_data.substr(0, colon_pos) : sample_data;
            
            GenotypeCall gt = parseGenotype(gt_str);
            
            if (!matchesCode(gt, sample.code)) {
                return false; // This variant doesn't match the genotype criteria
            }
            
            // Extract depth if depth filtering is enabled
            if (min_depth >= 0 || max_depth >= 0) {
                int depth = extractDepth(sample_data, format_field);
                if (depth >= 0) {
                    sample_depths.push_back(depth);
                }
            }
        }
        
        // Apply depth filtering if enabled
        if ((min_depth >= 0 || max_depth >= 0) && !sample_depths.empty()) {
            auto min_max = std::minmax_element(sample_depths.begin(), sample_depths.end());
            int min_sample_depth = *min_max.first;
            int max_sample_depth = *min_max.second;
            
            if (min_depth >= 0 && min_sample_depth < min_depth) {
                return false; // Minimum depth not met
            }
            
            if (max_depth >= 0 && max_sample_depth > max_depth) {
                return false; // Maximum depth exceeded
            }
        }
        
        return true; // All samples match their criteria
    }
};

std::string readLine(gzFile file) {
    char buffer[10000];
    if (gzgets(file, buffer, sizeof(buffer)) != nullptr) {
        std::string line(buffer);
        // Remove trailing newline
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
        }
        return line;
    }
    return "";
}

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options] -i <vcf_file> -s <sample_file>\n";
    std::cout << "Options:\n";
    std::cout << "  -i <file>    Input VCF file (.vcf or .vcf.gz)\n";
    std::cout << "  -s <file>    Sample file (format: sample_name code)\n";
    std::cout << "  -a           Allow missing data as neutral\n";
    std::cout << "  -m           Accept multiple alleles\n";
    std::cout << "  --min-depth <int>    Minimum depth filter (applied to all samples)\n";
    std::cout << "  --max-depth <int>    Maximum depth filter (applied to all samples)\n";
    std::cout << "  -o <file>    Output file (default: stdout)\n";
    std::cout << "  -h           Show this help message\n";
    std::cout << "\nSample codes:\n";
    std::cout << "  1 = Homozygous reference (0/0)\n";
    std::cout << "  2 = Homozygous alternative (1/1)\n";
    std::cout << "  3 = Heterozygous (0/1)\n";
    std::cout << "  4 = Heterozygous or homozygous reference\n";
    std::cout << "  5 = Homozygous alternative or heterozygous\n";
    std::cout << "  6 = Homozygous reference or homozygous alternative\n";
    std::cout << "  7 = Any genotype\n";
    std::cout << "\nDepth filtering:\n";
    std::cout << "  Depth filters are applied to all samples in the sample list.\n";
    std::cout << "  --min-depth: Minimum depth - the sample with lowest depth must have at least this coverage\n";
    std::cout << "  --max-depth: Maximum depth - the sample with highest depth must have at most this coverage\n";
}

std::vector<Sample> readSampleFile(const std::string& filename) {
    std::vector<Sample> samples;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open sample file " << filename << std::endl;
        return samples;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string sample_name;
        int code;
        
        if (iss >> sample_name >> code) {
            if (code >= 1 && code <= 7) {
                samples.push_back({sample_name, code});
            } else {
                std::cerr << "Warning: Invalid code " << code << " for sample " << sample_name << std::endl;
            }
        }
    }
    
    return samples;
}

int main(int argc, char* argv[]) {
    std::string vcf_file;
    std::string sample_file;
    std::string output_file;
    bool allow_missing = false;
    bool multiple_alleles = false;
    int min_depth = -1;
    int max_depth = -1;
    
    // Define long options for depth filtering
    static struct option long_options[] = {
        {"min-depth", required_argument, 0, 1},
        {"max-depth", required_argument, 0, 2},
        {0, 0, 0, 0}
    };
    
    int opt;
    int option_index = 0;
    
    while ((opt = getopt_long(argc, argv, "i:s:o:amh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i':
                vcf_file = optarg;
                break;
            case 's':
                sample_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'a':
                allow_missing = true;
                break;
            case 'm':
                multiple_alleles = true;
                break;
            case 'h':
                printUsage(argv[0]);
                return 0;
            case 1: // --min-depth
                min_depth = std::atoi(optarg);
                break;
            case 2: // --max-depth
                max_depth = std::atoi(optarg);
                break;
            default:
                printUsage(argv[0]);
                return 1;
        }
    }
    
    if (vcf_file.empty() || sample_file.empty()) {
        std::cerr << "Error: Both VCF file (-i) and sample file (-s) are required\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Read sample information
    std::vector<Sample> samples = readSampleFile(sample_file);
    if (samples.empty()) {
        std::cerr << "Error: No valid samples found in sample file\n";
        return 1;
    }
    
    // Initialize filter
    VCFFilter filter(samples, allow_missing, multiple_alleles, min_depth, max_depth);
    
    // Open VCF file (handles both .vcf and .vcf.gz)
    gzFile vcf = gzopen(vcf_file.c_str(), "r");
    if (!vcf) {
        std::cerr << "Error: Could not open VCF file " << vcf_file << std::endl;
        return 1;
    }
    
    // Open output file or use stdout
    std::ostream* output = &std::cout;
    std::ofstream outfile;
    if (!output_file.empty()) {
        outfile.open(output_file);
        if (!outfile.is_open()) {
            std::cerr << "Error: Could not open output file " << output_file << std::endl;
            gzclose(vcf);
            return 1;
        }
        output = &outfile;
    }
    
    std::string line;
    bool header_processed = false;
    
    // Process VCF file
    while (!(line = readLine(vcf)).empty()) {
        if (line.empty()) continue;
        
        if (line[0] == '#') {
            // Header line
            *output << line << std::endl;
            if (line.substr(0, 6) == "#CHROM") {
                filter.setSamplePositions(line);
                header_processed = true;
            }
        } else {
            // Data line
            if (!header_processed) {
                std::cerr << "Error: No header line found in VCF file\n";
                break;
            }
            
            if (filter.filterVariant(line)) {
                *output << line << std::endl;
            }
        }
    }
    
    gzclose(vcf);
    if (outfile.is_open()) {
        outfile.close();
    }
    
    return 0;
}
