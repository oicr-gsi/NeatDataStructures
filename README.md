# NeatDataStructures
This is the repository for code which uses human genomic data 
(specifically, genomic variants pooled from many individuals) to create data structures for the NEAT genomic read simulator.



## Instructions to grab human genomic data from ICGC
Access the ICGC Data Portal at https://dcc.icgc.org/. Select "Cancer Projects" and choose the cancer of interest. Data are organized by
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



## FindNucleotideContextOnReference.pl 

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
Each tumor_genotype is split and compared to the reference allele from the references genome. If the original germline allele does not match the reference genome allele, the middle allele of the trinucleotide context is replaced with the original germline allele.

### Generating Trinucleotide Frequency Matrices
The script creates a hash of hashes as it runs through the input cancer data. A hash is created for each trinucleotide context, when the middle allele is removed (N_N). Each trinucleotide hash contains hashes of the original mutated from allele, which each contain keys of the mutated to alleles. As the script runs through the input file, each key is tallied up.

trinucleotide_context{mutated_from}{mutated_to}

After running through the input file and tallying up the total of each key, in reference to the mutated from allele and trinucleotide context, each value is divided by the total number of it corresponding mutated from allele in the trinucleotide context. This produces a new matrix for each trinucleotide context, where the probabilities of each row add up to 1. Each matrix file is named "Breast_N-N.trinuc" corresponding to the trinucleotide context.

### Generating Insertion and Deletion Length Frequencies
As the script is running through the input file, it locates insertions and deletions and returns their lengths using a hash. For example, the insertion length hash has lengths of the insertions as its key
and the total number of each length as the values.

Insertion length totals are then divided by the total number of insertions to calculate the probabity of each insertion length occurring in the data set. These probabilities are written to the files "_insLength.prob" and "_delLength.prob".

### Generating Overall Mutation Type Frequencies
The script tallies up the number of insertions and deletions found in the input cancer data file, as well as the total number of mutations. Insertion and deletion totals are divided by the total number of mutations to calculate their probability in the data set. More mutation type probabilities are to be included later. These probabilities are written to the file "Breast_overall.prob".

### Generating Consequence Type Frequencies
ICGC tsv files contain a "Consequence Type" column that contains the consequence of each mutation. As the script runs through the input file, each unique consequence type is tallied up, and each total is divided by the total mutations in the data set.

## FindNucleotideContextOnReference.healthy.pl

###Overview
This script was created to extract the same frequency models from dbSNP vcfs as were extracted from ICGC tsv files. dbSNP is an NCBI database of healthy variants that is open for use and analysis by the public. NEAT is able to use non-cancer variants for simulation.

### Changes to InDel Recognition
For all InDels in dbSNP, the first nucleotide upstream of the insertion or deletion is included in both the REF(mutated from) and ALT(mutated to) columns. This required an edit to the script in the way it recognizes InDels. The healthy variant script defines the length of the REF and ALT columns in nucleotides, and recognizes when length(REF) does not equal length(ALT). Mutation length totals are calculated using an insertion hash and a deletion hash with length as the keys and total as the values.

### Analyzing Multi-Nucleotide Variants (in progress)
dbSNP contains a small, but significant, number of nucleotide substitutions greater than one base in length. These variants were not found in ICGC simple somatic mutation data, therefore an edit is required for the healthy variant script. This is soon to be implemented.

