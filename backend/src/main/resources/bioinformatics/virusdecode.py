from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from urllib.error import HTTPError
import subprocess
import os
import sys
import json
from io import StringIO
import requests
Entrez.email = "your_email@example.com"

current_dir = os.path.dirname(os.path.abspath(__file__))
# JSON 데이터를 파일로 저장하는 함수
def save_json(data, file_path):
    file_path = current_dir +"/data/"+ file_path
    with open(file_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)

# JSON 파일로부터 특정 값을 가져오는 함수
def get_json(file_path):
    file_path = current_dir +"/data/"+  file_path
    with open(file_path, 'r') as json_file:
        data = json.load(json_file)
        return data

def get_metadata(reference_id):
    # Get reference sequence
    try:
        handle = Entrez.efetch(db="nucleotide", id=reference_id, rettype="gb", retmode="text")
        reference_sequence = SeqIO.read(handle, "genbank")
        handle.close()
        metadata={
            "Sequence ID": reference_sequence.id,
            "Name": reference_sequence.name,
            "Description": reference_sequence.description,
            "Length": len(reference_sequence)
        }
        return metadata
    except HTTPError as e:
        print(f"HTTPError: {e.code} - {e.reason}")

class SequenceAlignment:
    def __init__(self, variant_sequences, reference_id, muscle_exe="muscle"):
        self.muscle_exe = muscle_exe
        self.reference_sequence = None
        self.variant_sequences=variant_sequences
        self.aligned_sequences = []
        self.reference_id = None
        self.alignment_dict={}
        self.protein_length={}
        self.reference_protein_seq=""
        self.alignment_index={}
        
        try:
            # Get reference sequence
            handle = Entrez.efetch(db="nucleotide", id=reference_id, rettype="gb", retmode="text")
            self.reference_sequence = SeqIO.read(handle, "genbank")
            handle.close()
            self.reference_id = self.reference_sequence.id

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
        # Run muscle using subprocess, reading from stdin and writing to stdout
        result = subprocess.run([self.muscle_exe], input=self.combined_memory_file.read(), text=True, capture_output=True)  # 수정된 부분: muscle 실행을 위해 stdin으로 데이터를 전달하고 stdout을 메모리에 저장
        if result.returncode == 0:
            # print("Muscle ran successfully")
            self.aligned_memory_file = StringIO(result.stdout)  # 수정된 부분: 결과를 StringIO 객체에 저장하여 메모리 내에서 처리
        # else:
            # print("Error running muscle:", result.stderr)

    def read_alignment(self):
        # Read alignment
        alignment = list(SeqIO.parse(self.aligned_memory_file, "fasta"))  # 수정된 부분: aligned_memory_file로부터 메모리 내 데이터를 읽어옴

        # Sort alignment
#         desired_order = [self.reference_sequence.id] + [record.id for record in self.variant_sequences]
        desired_order = [self.reference_sequence.id] + list(self.variant_sequences.keys())
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

    def get_alignment_index(self):
        return self.alignment_index

    def get_aligned_sequences(self):
        aligned_sequences_dict = {record.id: str(record.seq) for record in self.aligned_sequences}
        return aligned_sequences_dict


    def run(self):
        self.read_sequences()
        self.run_muscle_dna()
        self.read_alignment()
    



class SequenceAnalysis:
    def __init__(self, alignment_index, alignment_dict, reference_id, gene, variant_id, start, end):
        self.alignment_index = alignment_index
        self.alignment_dict = alignment_dict
        self.reference_id = reference_id
        self.gene = gene
        self.variant_id = variant_id
        self.start = start
        self.end = end
        self.target_sequence = None
        self.linearDesign = []
        self.protParam = []

    def run_linear_design(self):
        start=self.start
        end=self.end

        # Update region
        (idx_start,idx_end) = self.alignment_index[self.gene]
        input_sequence = self.alignment_dict[self.reference_id][idx_start:idx_end]

        gap_count = input_sequence[:start].count("-")
        start += gap_count
        end += gap_count

        gap_count = input_sequence[start:end].count("-")
        end += gap_count

        # Get target sequence
        input_sequence = self.alignment_dict[self.variant_id][idx_start:idx_end]
        input_sequence = input_sequence[start:end].replace("-", "")
        self.target_sequence = input_sequence
        
        # Run LinearDesign

        # Execute the command and capture the result
        os.chdir(os.path.join(current_dir, "LinearDesign"))
        command = f"echo {input_sequence} | ./lineardesign"
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        os.chdir(current_dir)

        # Check if the command was executed successfully
        if process.returncode == 0:
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

        # else:
            # print("Error executing command")
            # print(stderr)

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

    def get_linearDesign(self):
        mRNA_sequence, mRNA_structure, free_energy, cai = self.linearDesign
        linearDesign_dict = {
            "mRNA_sequence": mRNA_sequence,
            "mRNA_structure": mRNA_structure,
            "free_energy": free_energy,
            "cai": cai
        }
        return linearDesign_dict

    def get_protParam(self):
        sequence, molecular_weight, amino_acid_count, amino_acid_percent, isoelectric_point, instability_index, secondary_structure_fraction, gravy, aromaticity = self.protParam
        protParam_dict = {
            "sequence": sequence,
            "molecular_weight": molecular_weight,
            "amino_acid_count": amino_acid_count,
            "amino_acid_percent": amino_acid_percent,
            "isoelectric_point": isoelectric_point,
            "instability_index": instability_index,
            "secondary_structure_fraction": secondary_structure_fraction,
            "gravy": gravy,
            "aromaticity": aromaticity
        }
        return protParam_dict

    def run(self):
        self.run_linear_design()
        self.set_protParam()




if __name__ == "__main__":
    option = int(sys.argv[1])

    # metadata
    if option == 1:
        reference_id = sys.argv[2]

        os.makedirs(current_dir+"/data", exist_ok=True)

        metadata = get_metadata(reference_id)
        
        save_json(metadata, "metadata.json")  # JSON 파일로 저장
        print(json.dumps(metadata))

    # alignment
    elif option == 2:
        variant_sequences = {}

        # read fasta file
        for record in SeqIO.parse(current_dir+"/data/combined.fasta", "fasta"):
            variant_sequences[record.id] = record.seq

        # get metadata
        metadata = get_json("metadata.json")
        reference_id = metadata.get("Sequence ID", None)

        # run alignment
        alignment = SequenceAlignment(variant_sequences, reference_id)
        alignment.run()

        # get alignment data
        alignment_index = alignment.get_alignment_index()
        aligned_sequences_dict = alignment.get_aligned_sequences()
        alignment_data = {
            "alignment_index": alignment_index,
            "aligned_sequences": aligned_sequences_dict,
        }

        save_json(alignment_data, "alignment_data.json")  # JSON 파일로 저장
        print(json.dumps(alignment_data))

    # linearDesign, protparam data
    elif option == 3:
        # get metadata
        metadata = get_json("metadata.json")
        reference_id = metadata.get("Sequence ID", None)
        
        # get alignment data
        alignment_data = get_json("alignment_data.json")
        alignment_index = alignment_data.get("alignment_index", None)
        alignment_dict = alignment_data.get("aligned_sequences", None)
        
        # set gene, variant_id, start, end
        gene=sys.argv[2]
        variant_id=sys.argv[3]
        start = int(sys.argv[4])
        end = int(sys.argv[5])

        # run sequence analysis
        analysis = SequenceAnalysis(alignment_index, alignment_dict, reference_id, gene, variant_id, start, end)
        analysis.run()

        # get linearDesign and protParam data
        linearDesign_dict = analysis.get_linearDesign()
        protParam_dict = analysis.get_protParam()
        linearDesign_data = {
            "linearDesign": linearDesign_dict,
            "protParam": protParam_dict
        }

        save_json(linearDesign_data, "linearDesign_data.json")  # JSON 파일로 저장
        print(json.dumps(linearDesign_data))


    elif option == 4:
        # RCSB PDB API endpoint
        url = "https://search.rcsb.org/rcsbsearch/v2/query"

        # get metadata
        metadata = get_json("metadata.json")
        reference_id = metadata.get("Sequence ID", None)

        query = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": reference_id
                }
            },
            "return_type": "entry"
        }

        # request API
        response = requests.post(url, json=query)

        # get PDB IDs (list)
        if response.status_code == 200:
            results = response.json()
            pdb_list = [entry['identifier'] for entry in results['result_set']]
            pdb_dict = {f"value{idx+1}": pdb for idx, pdb in enumerate(pdb_list)}
            print(json.dumps(pdb_dict, indent=4))
        # else:
        #     print(f"Error: {response.status_code}")



