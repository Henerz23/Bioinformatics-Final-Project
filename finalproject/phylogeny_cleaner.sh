#!/bin/bash

##################################################
### Format sequence file for MSA and phylogeny ###
##################################################

# use awk to print the reformat the file to display only the species name
# this leaves only relevent information for phylogenetic analysis
awk '
/^>/ {
    print ">" $2 "_" $3;
    next;
}
{ print }
' translated_sequences.fasta > tmp.fasta

# remove any X amino acids to allow for cleaner MSA analysis 
sed '/^>/! s/X//g' tmp.fasta > translated_sequences.fasta
# remove unecessary file
rm tmp.fasta

#update user 
echo "Your MSA and phylogeny ready sequences are complete!"