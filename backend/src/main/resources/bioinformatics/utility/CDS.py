from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Entrez

Entrez.email = "your_email@example.com"

# Define the input and output file paths
output_fasta = "translated_protein_sequences.fasta"

# Define the region to translate (266 to 21555)
reference_id = "NC_045512.2"

# Read the DNA sequences from the input FASTA file
handle = Entrez.efetch(db="nucleotide", id=reference_id, rettype="gb", retmode="text")
dna_sequences = SeqIO.read(handle, "genbank")
handle.close()
# Create an empty list to hold the translated protein sequences
protein_sequences = []



locations = []
for feature in dna_sequences.features:
    if(feature.type == "CDS"):
        for part in feature.location.parts:
            start = part.start
            end = part.end
            locations.append((int(start), int(end),feature.qualifiers["gene"][0]))


# print(dna_sequences.annotations)
for annotation in dna_sequences.annotations:
    print(annotation, dna_sequences.annotations[annotation])
# print(locations)
# for location in locations:
#     print(location)

for start_position, end_position, gene in locations:
    region_seq = dna_sequences.seq[start_position:end_position]
    protein_seq = region_seq.translate(to_stop=True)
    protein_record = SeqRecord(protein_seq, id=gene, description=f"translated protein {start_position}-{end_position}")
    protein_sequences.append(protein_record)


SeqIO.write(protein_sequences, output_fasta, "fasta")

print(f"Translated protein sequences have been saved to {output_fasta}")
