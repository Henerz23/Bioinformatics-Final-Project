##################################################
####         Translate longest ORF's           ###
##################################################

# Import modules
# import glob and chain modules to loop over files in a directory
import glob
from itertools import chain
# import Bio python modules for efficient translation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq



## PART 1 ##
# Loop over the files in the directory to process samples as well as homolog sequences
sequences = []
# loop over .txt files downloaded from blast and .fasta files generated from samples in the merged sequences directory
patterns = ["merged_sequences/*.txt", "merged_sequences/*.fasta"]
for filepath in chain.from_iterable(glob.glob(p) for p in patterns):
    # add each sequence to list of sequences
    sequences.extend(list(SeqIO.parse(filepath, "fasta")))
    # update user
    print("Processing " + filepath)



## PART 2 ##
# Translate the sample and homolog sequences in 3 ORF's , and save the longest
# Store output amino-acid SeqRecord objects
aa_sequences = []  

# loop over the newly generated file of sequences
for record in sequences:
    
    # Nucleotide sequence from the FASTA record
    dna_seq = record.seq   

    # generate variable to store the longest ORF
    best_orf = ""
    
    # generate variable to store which frame contains the longest ORF
    best_frame = None

    # Translate sequence in all three forward reading frames
    for frame in range(3):
        
        # Translate starting at this frame; keep stop codons
        # utilise table 2, as the sequences are of vertabrate mitochondrial DNA
        protein = dna_seq[frame:].translate(table=2, to_stop=False)

        # Split translation into fragments separated by stop codons ("*")
        # Each fragment is a potential ORF (no internal stops)
        fragments = str(protein).split("*")

        # Find the longest continuous amino-acid stretch (longest ORF) in this frame
        longest_in_frame = max(fragments, key=len)

        # If this ORF is longer than the best one we've seen so far, store it
        # this will generate the most accurate MSA for species identification
        if len(longest_in_frame) > len(best_orf):
            best_orf = longest_in_frame
            best_frame = frame

    # Convert the longest ORF string back into a Seq object
    best_protein = Seq(best_orf)

    # Extract species name (Latin binomial) from the FASTA description
    # Assumes format: >ID Genus species ... 
    start = record.description.find("[organism=") + len("[organism=")
    end = record.description.find("]", start)
    latin_name = record.description[start:end].replace(" ", "_")  # replace space with underscore

    # Create a new SeqRecord containing the longest ORF for this sequence
    aa_record = SeqRecord(
        best_protein,
        id=record.id,                      # Keep the same FASTA ID
        name=latin_name,                   # Store species name
        description=record.description,    # Preserve full original header
        annotations={
            "type": "longest ORF",
            "frame": best_frame            # Which reading frame the ORF came from
        }
    )
    
    # update the user on the translation progress and the reading frame used
    print(aa_record.id + "'s sequence was translated in reading frame: " + str(best_frame + 1))
    
    # Add to output list
    aa_sequences.append(aa_record)

# Write all translated longest-ORF sequences to a FASTA file for efficient MSA
SeqIO.write(aa_sequences, "translated_sequences.fasta", "fasta")

