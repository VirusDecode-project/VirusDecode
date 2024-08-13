import os
import subprocess
from Bio import AlignIO, SeqIO

class SequenceAlignment:
    def __init__(self, file1, file2, muscle_exe="muscle"):
        self.file1 = file1
        self.file2 = file2
        self.muscle_exe = muscle_exe
        self.combined_file = "combined.fasta"
        self.aligned_file = "aligned.fasta"
        self.base_colors = {
            'A': '#ff9999',
            'T': '#99ff99',
            'G': '#9999ff',
            'C': '#ffff99',
            '-': '#cccccc'
        }
        self.seq1 = None
        self.seq2 = None
        self.seq1_aligned = None
        self.seq2_aligned = None
        self.mutations = []

    def read_sequences(self):
        self.seq1 = SeqIO.read(self.file1, "fasta")
        self.seq2 = SeqIO.read(self.file2, "fasta")

    def combine_sequences(self):
        with open(self.combined_file, "w") as f:
            SeqIO.write([self.seq1, self.seq2], f, "fasta")

    def run_muscle(self):
        subprocess.run([self.muscle_exe, "-in", self.combined_file, "-out", self.aligned_file])

    def read_alignment(self):
        alignment = AlignIO.read(self.aligned_file, "fasta")
        self.seq1_aligned = alignment[0].seq
        self.seq2_aligned = alignment[1].seq

    def find_mutations(self):
        self.mutations = [
            i for i, (a, b) in enumerate(zip(self.seq1_aligned, self.seq2_aligned))
            if a != b and b != 'N'
        ]
        # 돌연변이 정보 출력
        print(f"Total number of mutations: {len(self.mutations)}")
        print("Positions of mutations:")
        for mut in self.mutations:
            print(f"Position {mut + 1}: {self.seq1_aligned[mut]} -> {self.seq2_aligned[mut]}")

    def generate_html(self):
        html_content = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Sequence Alignment Mutations</title>
    <style>
        .base { display: inline-block; padding: 2px; margin: 1px; }
        .mutation { background-color: red; color: white; }
    </style>
</head>
<body>
    <h1>Sequence Alignment Mutations</h1>
    <pre>
"""
        for i, (a, b) in enumerate(zip(self.seq1_aligned, self.seq2_aligned)):
            if i % 60 == 0 and i != 0:
                html_content += "\n"
            if i in self.mutations:
                html_content += f'<span class="base mutation">{b}</span>'
            elif b == 'N':
                html_content += f'<span class="base" style="background-color: {self.base_colors["-"]};">{b}</span>'
            else:
                html_content += f'<span class="base" style="background-color: {self.base_colors[b]};">{b}</span>'
        
        html_content += """
    </pre>
</body>
</html>
"""
        with open("mutations.html", "w") as f:
            f.write(html_content)
        print("HTML file 'mutations.html' has been created.")

    def save_mutation_info(self):
        with open("mutation.txt", "w") as f:
            f.write(f"Total number of mutations: {len(self.mutations)}\n")
            f.write("Positions of mutations:\n")
            for mut in self.mutations:
                f.write(f"Position {mut + 1}: {self.seq1_aligned[mut]} -> {self.seq2_aligned[mut]}\n")
        print("Text file 'mutation.txt' has been created.")

    def run(self):
        self.read_sequences()
        self.combine_sequences()
        self.run_muscle()
        self.read_alignment()
        self.find_mutations()
        self.generate_html()
        self.save_mutation_info()

if __name__ == "__main__":
    alignment = SequenceAlignment("NC_045512.fasta", "OX014251.fasta")
    alignment.run()
