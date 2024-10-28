import re
import sys
import json
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import subprocess
Entrez.email = "your_email@example.com"

class Alignment:
    def __init__(self, reference_id, fasta_content, muscle_exe="muscle"):
        self.muscle_exe = muscle_exe
        self.reference_sequence = None
        self.variant_sequences={}
        self.fasta_content=fasta_content
        self.aligned_sequences = []
        self.reference_id = reference_id
        self.alignment_dict={}
        self.protein_length={}
        self.reference_protein_seq=""
        self.alignment_index={}
        self.alignment_data={}


    # # ATCG 문자만 허용하는 정규식 검사 함수
    def validate_sequences(self):
        def is_valid_sequence(sequence: str) -> bool:
            return bool(re.fullmatch(r'[ATCGatcg]*', str(sequence)))

        variant_sequences = {}
        fasta_io = StringIO(self.fasta_content)
        try:
            for record in SeqIO.parse(fasta_io, "fasta"):
                if not is_valid_sequence(record.seq):
                    sys.stderr.write(f"Invalid sequence for ${record.id}.\nOnly A, T, C, and G are allowed.")
                    sys.exit(22)
                variant_sequences[record.id] = record.seq
        except Exception as e:
            sys.stderr.write(f"fasta format error in SeqIO.parse.\n{str(e)}")
            sys.exit(22)
        self.variant_sequences = variant_sequences

    def fetch_reference_sequence(self):
        try:
            # Get reference sequence
            handle = Entrez.efetch(db="nucleotide", id=self.reference_id, rettype="gb", retmode="text")
            self.reference_sequence = SeqIO.read(handle, "genbank")
            handle.close()
            self.reference_id = self.reference_sequence.id

        except Exception as e:
            sys.stderr.write(f"There was an error fetching the reference sequence from NCBI.\n{str(e)}")
            sys.exit(11)

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
        translated_variants = {id: seq.translate(to_stop=True)
                               for id, seq in self.variant_sequences.items()}
        translated_records = [SeqRecord(seq=translation, id=id, description="translated variant protein")
                              for id, translation in translated_variants.items()]


        # Combine reference and translated variant sequences
        reference_protein_record = SeqRecord(Seq(self.reference_protein_seq), id=self.reference_sequence.id, description="reference protein")
        self.combined_memory_file = StringIO() 
        SeqIO.write([reference_protein_record] + translated_records, self.combined_memory_file, "fasta")
        self.combined_memory_file.seek(0)

    def run_muscle_dna(self):
        try:
            # Run muscle using subprocess, reading from stdin and writing to stdout
            result = subprocess.run(
                [self.muscle_exe], 
                input=self.combined_memory_file.read(), 
                text=True, 
                capture_output=True,
                check=True  # CalledProcessError를 발생시키도록 설정
            )
            self.aligned_memory_file = StringIO(result.stdout)  # 결과를 메모리에 저장하여 처리
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"Error running MUSCLE: {e}\n")
            sys.exit(21)

    def read_alignment(self):
        # Read alignment
        alignment = list(SeqIO.parse(self.aligned_memory_file, "fasta"))  # 수정된 부분: aligned_memory_file로부터 메모리 내 데이터를 읽어옴

        # Sort alignment
        desired_order = [self.reference_sequence.id] + list(self.variant_sequences.keys())
        self.alignment_dict = {record.id: record for record in alignment}
        length_of_first_seq = len(self.alignment_dict[self.reference_sequence.id].seq)
        for id in desired_order:
            if id not in self.alignment_dict:
                self.alignment_dict[id] = SeqRecord(Seq("-" * length_of_first_seq), id=id)
                
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

    def get_alignment_index(self):
        return self.alignment_index

    def get_aligned_sequences(self):
        aligned_sequences_dict = {record.id: str(record.seq) for record in self.aligned_sequences}
        return aligned_sequences_dict

    def set_alignment_data(self):
        alignment_index = self.get_alignment_index()
        aligned_sequences_dict = self.get_aligned_sequences()
        alignment_data = {
            "alignment_index": alignment_index,
            "aligned_sequences": aligned_sequences_dict,
        }
        self.alignment_data = alignment_data

    def print_alignment_data(self):
        print(json.dumps(self.alignment_data, indent=4))

    def run(self):
        self.validate_sequences()
        self.fetch_reference_sequence()
        self.read_sequences()
        self.run_muscle_dna()
        self.read_alignment()
        self.set_alignment_data()
        self.print_alignment_data()