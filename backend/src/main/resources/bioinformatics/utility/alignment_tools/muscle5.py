import subprocess
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class SequenceAlignment:
    def __init__(self, files, reference_file, muscle_exe="/usr/local/bin/muscle5"):
        self.reference_file = reference_file
        self.files = files
        self.files_len = len(files)
        self.muscle_exe = muscle_exe
        self.combined_file = "combined.fasta"
        self.aligned_file = "aligned.fasta"
        self.protein_file = "proteins.fasta"
        self.aligned_protein_file = "aligned_proteins.fasta"
        self.reference_sequence = None
        self.sequences = []
        self.aligned_sequences = []
        self.protein_sequences = []
        self.mutations = []

    def read_sequences(self):
        self.reference_sequence = SeqIO.read(self.reference_file, "fasta")
        self.sequences = [SeqIO.read(file, "fasta") for file in self.files]

    def combine_sequences(self):
        with open(self.combined_file, "w") as f:
            SeqIO.write([self.reference_sequence] + self.sequences, f, "fasta")
        # print("Combined file content:")
        # with open(self.combined_file, "r") as f:
            # print(f.read())

    def run_muscle_dna(self):
        result = subprocess.run([self.muscle_exe, "-align", self.combined_file, "-output", self.aligned_file])
        if result.returncode != 0:
            print("Error running MUSCLE:", result.stderr)
        else:
            print("MUSCLE alignment complete")
        print("Aligned file content:")
        with open(self.aligned_file, "r") as f:
            print(f.read())

    def read_alignment(self):
        alignment = AlignIO.read(self.aligned_file, "fasta")
        self.aligned_sequences = [record for record in alignment]
        print("len of aligned_seq", len(self.aligned_sequences))

    def translate_sequences(self):
        self.protein_sequences = []
        for record in self.aligned_sequences:
            clean_dna_seq = str(record.seq).replace('-', 'N')
            clean_dna_seq = Seq(clean_dna_seq)
            protein_seq = clean_dna_seq.translate()
            self.protein_sequences.append(SeqRecord(protein_seq, id=record.id))

    def write_protein_sequences(self):
        with open(self.protein_file, "w") as f:
            SeqIO.write(self.protein_sequences, f, "fasta")

    def find_mutation(self):
        print(self.files_len)
        for variant_id in range(self.files_len):
            mutation = [
                i for i, (a, b) in enumerate(zip(self.protein_sequences[0].seq, self.protein_sequences[variant_id + 1].seq))
                if a != b and b != 'X'
            ]
            self.mutations.append(mutation)
            print(f"Total mutations: {len(mutation)}/{len(self.protein_sequences[0].seq)}")

    def print_mutation(self):
        for variant_id in range(self.files_len):
            print(f"Mutations for variant {variant_id + 1}:")
            for mutation in self.mutations:
                for mut in mutation:
                    print(f"Position {mut + 1}: {self.protein_sequences[0].seq[mut]} -> {self.protein_sequences[variant_id + 1].seq[mut]}")

    def run(self):
        self.read_sequences()
        self.combine_sequences()
        self.run_muscle_dna()
        self.read_alignment()
        self.translate_sequences()
        self.write_protein_sequences()
        self.find_mutation()
        self.print_mutation()

if __name__ == "__main__":
    reference_file = "reference.fasta"
    files = ["variant2.fasta","variant3.fasta"]
    alignment = SequenceAlignment(files, reference_file)
    alignment.run()
