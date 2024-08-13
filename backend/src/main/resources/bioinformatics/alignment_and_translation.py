from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from urllib.error import HTTPError
import subprocess
import os
import sys
import json

current_dir = os.path.dirname(os.path.abspath(__file__))

class SequenceAlignment:
    def __init__(self, files, reference_id, muscle_exe="muscle"):
        Entrez.email = "your_email@example.com"
        self.files = files
        self.muscle_exe = muscle_exe
        self.combined_file = os.path.join(current_dir, "result/combined.fasta")
        self.aligned_file = os.path.join(current_dir, "result/aligned.fasta")
        self.reference_sequence = None
        self.variant_sequences = []
        self.aligned_sequences = []
        self.reference_id = None
        self.alignment_dict={}
        self.protein_length={}
        self.mutation_dict={}
        self.reference_protein_seq=""
        self.alignment_index={}
        self.target_sequence = None
        self.protParam = []
        self.linearDesign = []
        
        try:
            # Get reference sequence
            handle = Entrez.efetch(db="nucleotide", id=reference_id, rettype="gb", retmode="text")
            self.reference_sequence = SeqIO.read(handle, "genbank")
            handle.close()
            self.reference_id = self.reference_sequence.id

            self.metadata={
                "Sequence ID": self.reference_sequence.id,
                "Name": self.reference_sequence.name,
                "Description": self.reference_sequence.description,
                "Length": len(self.reference_sequence)
            }
        except HTTPError as e:
            print(f"HTTPError: {e.code} - {e.reason}")
            

    def read_sequences(self):
        # Get reference protein sequence
        CDS_dict={}
        for feature in self.reference_sequence.features:
            if feature.type == 'CDS':
                gene = feature.qualifiers["gene"][0]
                if gene not in CDS_dict:
                    CDS_dict[gene] = feature.qualifiers["translation"][0]
                
        for gene, protein_seq in CDS_dict.items():
            self.protein_length[gene] = len(protein_seq)
            self.reference_protein_seq += protein_seq

        # Get variant sequences
        self.variant_sequences = [SeqIO.read(file, "fasta") for file in self.files]
        translated_variants = [variant.seq.translate(to_stop=True) for variant in self.variant_sequences]
        translated_records = [SeqRecord(Seq(translation), id=variant.id, description="translated variant protein")
                              for variant, translation in zip(self.variant_sequences, translated_variants)]
        
        # Combine reference and translated variant sequences
        reference_protein_record = SeqRecord(Seq(self.reference_protein_seq), id=self.reference_sequence.id, description="reference protein")
        with open(self.combined_file, "w") as f:
            SeqIO.write([reference_protein_record] + translated_records, f, "fasta")

    def run_muscle_dna(self):
        # Run muscle
        result = subprocess.run(
            [self.muscle_exe, "-in", self.combined_file, "-out", self.aligned_file],
            stdout=subprocess.DEVNULL,  # 표준 출력을 숨김
            stderr=subprocess.DEVNULL   # 표준 오류를 숨김
        )
        # if result.returncode != 0:
        #     print("Error running muscle:")
        # else:
        #     print("Muscle ran successfully:")


    def read_alignment(self):
        # Read alignment
        alignment = list(SeqIO.parse(self.aligned_file, "fasta"))

        # Sort alignment
        desired_order = [self.reference_sequence.id] + [record.id for record in self.variant_sequences]
        self.alignment_dict = {record.id: record for record in alignment}
        self.aligned_sequences = [self.alignment_dict[id] for id in desired_order]

        # Update protein length
        record = self.aligned_sequences[0]
        start=0
        for gene, length in self.protein_length.items():
            end = start + length
            gap_count = record.seq[start:end].count('-')
            end += gap_count
            self.protein_length[gene] = length+gap_count

            # Set alignment data
            self.alignment_index[gene] = (start, end)
            start=end


    def write_protein_sequences(self):
        # for (gene, start, end) in self.alignment_index:
        for gene, (start, end) in self.alignment_index.items():
            protein_sequences=[]
            for record in self.aligned_sequences:
                protein_record = SeqRecord(record.seq[start:end], id=record.id, description=f"translated protein {start}-{end}")
                protein_sequences.append(protein_record)

            with open(os.path.join(current_dir, f"result/proteins_{gene}.fasta"), "w") as f:
                SeqIO.write(protein_sequences, f, "fasta")


    def set_mutation(self):
        reference_protein = self.aligned_sequences[0].seq
        for record in self.aligned_sequences[1:]:
            variant_protein = record.seq
            mutation = []
            for i, (ref, var) in enumerate(zip(reference_protein, variant_protein)):
                if ref != var and ref != "-" and var != "-":
                    mutation.append((i, ref, var))
            self.mutation_dict[record.id] = mutation

    def run_linear_design(self, gene, variant_id):
        # Set RBD region
        RBD_start = 318
        RBD_end = 541

        # Update RBD region
        (start,end) = self.alignment_index[gene]
        input_sequence = str(self.alignment_dict[self.reference_id].seq[start:end])

        gap_count = input_sequence[:RBD_start].count("-")
        RBD_start += gap_count
        RBD_end += gap_count

        gap_count = input_sequence[RBD_start:RBD_end].count("-")
        RBD_end += gap_count

        # Get RBD region
        input_sequence = str(self.alignment_dict[variant_id].seq[start:end])
        input_sequence = input_sequence[RBD_start:RBD_end].replace("-", "")
        self.target_sequence = input_sequence
        
        # Run LinearDesign
        os.chdir("./LinearDesign")
#         os.path.join(current_dir, "./LinearDesign")

        # Execute the command and capture the result
        command = f"echo {input_sequence} | ./lineardesign"
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        # Check if the command was executed successfully
        if process.returncode == 0:
            print("Command executed successfully")

            # Save the output result as a list of lines
            output_lines = stdout.decode().splitlines()
            mRNA_sequence = output_lines[-4].replace('mRNA sequence:', '').strip()
            mRNA_structure = output_lines[-3].replace('mRNA structure:', '').strip()
            parts = output_lines[-2].split(';')
            free_energy = parts[0].replace('mRNA folding free energy:', '').strip()
            cai = parts[1].replace('mRNA CAI:', '').strip()

            # Set the linear design data
            self.linearDesign.append(mRNA_sequence)
            self.linearDesign.append(mRNA_structure)
            self.linearDesign.append(free_energy)
            self.linearDesign.append(cai)

        else:
            print("Error executing command")
            
    def set_protParam(self):
        # Protein sequence to analyze
        sequence = self.target_sequence

        # Create a protein analysis object
        protein_analysis = ProteinAnalysis(sequence)

        # Calculate molecular weight
        molecular_weight = protein_analysis.molecular_weight()

        # Count of amino acids
        amino_acid_count = protein_analysis.count_amino_acids()

        # Amino acid percentage
        amino_acid_percent = protein_analysis.get_amino_acids_percent()

        # Isoelectric point (pI)
        isoelectric_point = protein_analysis.isoelectric_point()

        # Instability index
        instability_index = protein_analysis.instability_index()

        # Fraction of polar, nonpolar, basic, and acidic amino acids
        secondary_structure_fraction = protein_analysis.secondary_structure_fraction()

        # Grand average of hydropathicity (GRAVY)
        gravy = protein_analysis.gravy()

        # Proportion of aromatic residues
        aromaticity = protein_analysis.aromaticity()

        # Return the protein parameters
        self.protParam.append(sequence)
        self.protParam.append(molecular_weight)
        self.protParam.append(amino_acid_count)
        self.protParam.append(amino_acid_percent)
        self.protParam.append(isoelectric_point)
        self.protParam.append(instability_index)
        self.protParam.append(secondary_structure_fraction)
        self.protParam.append(gravy)
        self.protParam.append(aromaticity)

    def get_metadata(self):
        return self.metadata
    
    def get_alignment_data(self):
        return self.alignment_index, self.aligned_sequences

    def get_mutation(self):
        return self.mutation_dict
    
    def get_linearDesign(self):
        return self.linearDesign

    def get_protParam(self):
        return self.protParam
        

    def run(self):
        self.read_sequences()
        self.run_muscle_dna()
        self.read_alignment()
        # self.write_protein_sequences()    # 부가기능(LinearDesign) 활용 시 주석 해제
        self.set_mutation()
#         self.run_linear_design("S", "MW642250.1")
#         self.set_protParam()


def string_to_json(input_string):
    # JSON 데이터 생성
    data = {
        "message": input_string
    }
    return data


if __name__ == "__main__":
    # reference_id = sys.argv[1]

    reference_id = "NC_045512"
#     files = ["./data/OL672836.1.spike.fasta", "./data/MW642250.1.spike.fasta", "./data/OM958567.1.spike.fasta"]

    files = [
        os.path.join(current_dir, "data/OL672836.1.spike.fasta"),
        os.path.join(current_dir, "data/MW642250.1.spike.fasta"),
        os.path.join(current_dir, "data/OM958567.1.spike.fasta")
    ]

    alignment = SequenceAlignment(files, reference_id)


    # alignment 실행 전 metadata 받아오기
    metadata = alignment.get_metadata()
#     for key, value in metadata.items():
#         print(f"{key}: {value}")
#     print(json.dumps(metadata))


    # alignment 실행
    alignment.run()
#
#     # alignment data 받아오기
    alignment_index, aligned_sequences = alignment.get_alignment_data()
    aligned_sequences_dict = {record.id: str(record.seq) for record in aligned_sequences}
    # for gene, (start, end) in alignment_index.items():
    #     print(gene, start, end)
    # for record in aligned_sequences:
    #     print(record.id)
    #     print(record.seq)
    #     print()

    # mutation data 받아오기
    mutation_dict = alignment.get_mutation()
    # for key, value in mutation_dict.items():
    #     print(key)
    #     for i, ref, var in value:
    #         print(f"{i}: {ref} -> {var}")
    #     print()

    alignment_data = {
        "alignment_index": alignment_index,
        "aligned_sequences": aligned_sequences_dict,
        "mutation_data": mutation_dict
    }
    print(json.dumps(alignment_data))

#
#     # linearDesign data 받아오기
#     linearDesign = alignment.get_linearDesign()
#     mRNA_sequence, mRNA_structure, free_energy, cai = linearDesign
#     linearDesign_dict = {
#         "mRNA_sequence": mRNA_sequence,
#         "mRNA_structure": mRNA_structure,
#         "free_energy": free_energy,
#         "cai": cai
#     }
#     # print(f"mRNA sequence: {mRNA_sequence}")
#     # print(f"mRNA structure: {mRNA_structure}")
#     # print(f"mRNA folding free energy: {free_energy}")
#     # print(f"mRNA CAI: {cai}")
#
#     # protparam data 받아오기
#     protParam = alignment.get_protParam()
#     sequence, molecular_weight, amino_acid_count, amino_acid_percent, isoelectric_point, instability_index, secondary_structure_fraction, gravy, aromaticity = protParam
#     protParam_dict = {
#         "sequence": sequence,
#         "molecular_weight": molecular_weight,
#         "amino_acid_count": amino_acid_count,
#         "amino_acid_percent": amino_acid_percent,
#         "isoelectric_point": isoelectric_point,
#         "instability_index": instability_index,
#         "secondary_structure_fraction": secondary_structure_fraction,
#         "gravy": gravy,
#         "aromaticity": aromaticity
#     }
#
#     # print(f"Protein Sequence: {sequence}")
#     # print(f"Molecular Weight: {molecular_weight:.2f} Da")
#     # print("Amino Acid Count:")
#     # for aa, count in amino_acid_count.items():
#     #     print(f"{aa}: {count}")
#     # print("Amino Acid Percent:")
#     # for aa, percent in amino_acid_percent.items():
#     #     print(f"{aa}: {percent:.2%}")
#     # print(f"Isoelectric Point (pI): {isoelectric_point:.2f}")
#     # print(f"Instability Index: {instability_index:.2f}")
#     # print(f"Secondary Structure Fraction (Helix, Turn, Sheet): {secondary_structure_fraction}")
#     # print(f"Gravy: {gravy:.2f}")
#     # print(f"Aromaticity: {aromaticity:.2%}")
#
#     linearDesign_data = {
#         "linearDesign": linearDesign_dict,
#         "protParam": protParam_dict
#     }
#     print(json.dumps(linearDesign_data))

