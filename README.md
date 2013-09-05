# CMSC737 - Assembly Model Testing #

## Goal ##
Design software to determine whether an assembly invalidates a set of constraints.

## Motivation ##
The genome sequence of an organism is a vital resource for biologists trying to better understand its function and evolution. Generating this sequence is not an easy task as modern sequencing technologies can only “read” small pieces of the genome. These sequences, known as reads, have to be pieced together by tools called assemblers using a collection of different heuristics since in almost all practical cases, **assemblers** cannot fully and accurately reconstruct the genome.

## Background ##
*READ FIRST!* http://www.cbcb.umd.edu/research/assembly_primer.shtml.  This will give you a better background than I ever could.

### Sequence/read ###
A single DNA sequence produced from some experiment.
http://en.wikipedia.org/wiki/DNA_sequencing

### FASTA ###
A flat file to store DNA sequences in the format of “>header [newline] dna sequence”.  *The majority of your data will be in this format!*
```
>header_name
ATCGTCAT
AGGATACA
>second_header
TTTCATTCC
```
http://en.wikipedia.org/wiki/FASTA_format

### Mate-pair ###
In most sequencing projects, the sizes of the fragments generated through the shotgun process are carefully controlled, thus providing a link between the sequence reads generated from the ends of a same fragment (called paired ends or mate pairs). 

### De novo assembly ###
Assembling the data by only looking at the reads without using any prior known information about the genome (opposite of referenced-based assembly).

### Contig ###
Contiguous pieces of DNA.
http://en.wikipedia.org/wiki/Contig

### SAM (Sequence Alignment/Map) format ###
SAM format is a generic format for storing sequence alignments, e.g., from reads to assemblies/contigs/references/etc.

http://samtools.sourceforge.net/

http://samtools.sourceforge.net/SAMv1.pdf

 
## Constraints ##

### Read coverage ###
#### Input ####
#### Output ####

### Mate-pair size agreement ###
#### Input ####
#### Output ####
