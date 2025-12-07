#!/bin/bash

##################################################
#### Make FASTQ samples into clean FASTA files ###
##################################################

# PART 1 #
# Use a for loop to clean the FASTQ files before formatting
# *requires prior identification of malformatted samples
# loop over files which are parts of the same sample using their distinct letters 
for letters in {A..D}; do

	# remove any previous files of the same name
	rm "sample${letters}.FASTQ"	
	
	# First, loop over the different parts of each letter to enable part specific control
	for f in sample${letters}_part*.FASTQ; do
		
		
		## FORMATTING TITLES ##
		# update user on processing progression
		echo "Cleaning: $f"

		# use awk to delete repeated lines in each file
		awk -i inplace '!seen[$0]++' "$f"
		
		# change any spaces for underscores in the FASTQ titles to make them the identical
		sed -i '1s/ /_/g' "$f"
		
		# delete empty lines
		sed -i '/^$/d' "$f"
		
		# replace titles with different structures to ensure they are universal
		sed -i 's/part_/part/g' "$f"

		# concatenation of all the same samples (E.g. A, B, C, D)
		# add a space between lines to ensure separate parts are concatenated on separate lines
		echo >> "sample${letters}.FASTQ"
		
		# add f into a new file called 
		cat "$f" >> "sample${letters}.FASTQ"
	done
	

	# PART 2 #
	# Convert the files from FASTQ to FASTA format, accounting for low Q scores
	# Then move the files into directories for further analysis
	
	
	# update the person on the progress
	echo "Converting Sample: $letters into FASTA format and concatenating"
	echo " "
	
	# create a file from the combined samples file for further editing
	file="sample${letters}.FASTQ"
	
	# delete empty lines
	sed -i '/^$/d' "$file"
	
	# Create one of two output FASTA files:
	# create a file where the parts are kept separate
	# this allows for analysis by parts of homologs, to check for anomilies
	separate="${file%.FASTQ}.fasta"
	

	# use seqtk from github to convert the files from FASTQ to FASTA in 3 stages:
	# 1- replace any bases with a Q score less than 99% into Ns
	# 2- convert the files to fasta format
	# 3- delete any spaces in the sequence
	seqtk seq -q20 -n N "$file" | seqtk seq -a | tr -d ' ' > "$separate"

	# create the second of 2 output files
	# this is a file where the parts are merged into one sequence
	# this enables analysis of homologues as a whole sequence
	merged="Merged_${file%.FASTQ}.fasta"
	
	# input a general header into the merged file
	echo ">sample${letters} Unknown species_${letters}" > "$merged"

	# delete any headers and spaces from the separate file and output this into the merged file
	sed '/^>/d; s/ //g' "$separate" >> "$merged"

	# move the merged file to the desired directory for further analysis
	mv "$merged" ../merged_sequences/
	
	# move the separate file into a distinct directory for analysis of parts separately
	# this aids with anomaly detection
	mv "$separate" ../separate_sequences
done

# Update user on file location and process completion
echo "Concatenated samples have been added to merged_sequences/"
echo "Separate samples have been added to separate_sequences/"
echo " "
echo "Cleaning and FASTA conversion complete!"
