# Variant-Filtering-by-GT-and-DP (VFGD)
A tool to filter a multi-sample VCF by genotype and depth of coverage.
____________________________________________________________________________________________________

Use at your own risk. I cannot provide support. All information obtained/inferred with this script is without any implied warranty of fitness for any purpose or use whatsoever.

ABOUT: 

This program takes a multi-sample VCF (compressed or not compressed) as input along with a two column sample list that contains a code of acceptable genotypes for each sample.  The VCF is then filtered for variants that meet the supplied criteria. Additional parameters can be provided for more refined filtering. The program reads the VCF header to map sample names to column positions. For each variant line, it checks if all specified samples match their filter code. Only outputs variants where all samples meet their criteria.

The command line options are:

-a: Allow missing data (./.) as neutral

-m: Accept multiple alleles (e.g., 2/2, 0/2, 1/3, etc.)

--min-depth: integer min depth

--max-depth: integer maximum depth 

-i: Input VCF file

-s: Sample file

-o: Output file (optional, defaults to stdout)


The sample file contains the sample name as it appears in the VCF in the first column and a numeric code in the second that defines the allowable genotype state:

   1 = Homozygous reference (0/0)
   
   2 = Homozygous alternative (1/1)
   
   3 = Heterozygous (0/1)
   
   4 = Heterozygous or homozygous reference
   
   5 = Homozygous alternative or heterozygous
   
   6 = Homozygous reference or homozygous alternative
   
   7 = Any genotype  

Filter Code Logic:

When -m is used, homozygous alternative accepts any identical non-zero alleles (1/1, 2/2, 3/3).
Heterozygous accepts any combination where one allele is 0.
Missing data (./.) is accepted as neutral only when -a flag is used.


PREREQUISITES:

1. C++ compiler (for example GCC/G++)
   
INSTALLATION: 

This program should work on all systems. Download the .cpp file and compile according to your system. For Linux Ubuntu, compile with g++ (g++ -o softmask_bed g++ -o vcf_filter VFGD_V1_5.cpp -lz). 

TO RUN:

Below are some usage examples:  


Basic filtering

./vcf_filter -i input.vcf.gz -s samples.txt

Allow missing data and multiple alleles

./vcf_filter -i input.vcf.gz -s samples.txt -a -m -o filtered.vcf

With depth

./vcf_filter -i test1.vcf.gz -s samples.txt -a -m --min-depth 40 --max-depth 200 -o filtered4depth.vcf

Show help

./vcf_filter -h
