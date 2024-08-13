from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Entrez

Entrez.email = "your_email@example.com"

# Define the region to translate (266 to 21555)
reference_id = "OM958567.1"

# Define the input and output file paths
output_fasta = f"{reference_id}.spike.fasta"


# Read the DNA sequences from the input FASTA file
handle = Entrez.efetch(db="nucleotide", id=reference_id, rettype="gb", retmode="text")
dna_sequences = SeqIO.read(handle, "genbank")
handle.close()

for feature in dna_sequences.features:
    if(feature.type == "CDS"):
        for part in feature.location.parts:
            start = part.start
            end = part.end
            gene = feature.qualifiers['gene'][0]
            if(feature.qualifiers['gene'][0] == "S"):
                # print(f"Gene: {gene}")
                # print(f"Start: {start}")
                # print(f"End: {end}")
                # print(f"Sequence: {dna_sequences.seq[start:end]}")
                # print()
                cutted_sequence = SeqRecord(dna_sequences.seq[start:end], id=dna_sequences.id, description=f"{gene} CDS:{start}-{end}")
                SeqIO.write(cutted_sequence, output_fasta, "fasta")
                break

