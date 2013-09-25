#! /usr/bin/bash

# First build the index that bowtie2 will use the align the reads.
bowtie2-build ../rhodobacter.assembly.fasta rhodobacter.bt2

# Run bowtie2 to get the mate pair alignments.  400 and 600 are the upper and lower bounds to search for a proper mate pair alignment.
bowtie2 -x rhodobacter.bt2 -f -1 ../rhodobacter.sequences.1.fasta -2 ../rhodobacter.sequences.2.fasta -I 400 -X 600 -S rhodobacter.sequences.sam