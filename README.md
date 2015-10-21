# NeatDataStructures
This is the repository for code which uses human genomic data 
(specifically, genomic variants pooled from many individuals) to create data structures for the NEAT genomic read simulator.



## Instructions to grab human genomic data from ICGC
Access the ICGC Data Portal at https://dcc.icgc.org/. Select "Cancer Projects" and choose the cancer of interest (in this case, Breast Cancer). Data are organized by
several categories that can be selected:
* Country
* Available Data Types
  * simple somatic mutations
  * copy number somatic mutations
  * structural somatic mutations
  * simple germline variants
  * array-based DNA-methylation
  * sequencing-based DNA-methylation
  * array-based gene expression
  * sequencing-based gene expression
  * protein expression
  * sequence-based miRNA expression
  * exon junctions
* Tumor Type

Once desired categories are selected, click on "Download Donor Data" to show file sizes and confirm file download.

## Controlled Data and Germline-Reference Allele Mismatch Information

For this work, comparisons of the variant statistics in tumor relative to the control tissue are desired.
Therefore, one needs to have access to the germline data, as follows.

ICGC's "Access Controlled Data" documention can be found at http://docs.icgc.org/access-controlled-data. To have access to controlled germline data, a DACO must be
submitted. Open tier data can be obtained without a DACO, but germline alleles that do not match the reference genome are masked and replaced with the reference
allele. Controlled data includes unmasked germline alleles.



## FindNucleotideContextOnReference 

###Overview
This script takes in human genomic data from the ICGC data portal as a TSV file. The coordinate for each variant in TSV file is located within 
the HG19 human reference. The corresponding trinucleotide context around that location on the reference is returned into a new column. 


### Running the Script
The script requires 4 arguments to be entered after the full path to FindNucleotideContextOnReference.pl

1. Full path to Fastahack
2. Full path to Reference Genome
3. Full path to input file
4. Full path to output file

### Computing the germline-tumor Allele Mismatch Information

TBC

