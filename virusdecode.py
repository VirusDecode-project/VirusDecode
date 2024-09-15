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
        print("NCBI에 요청한 nucleotide ID가 존재하지 않습니다.")
        sys.exit(11)

def get_pdb_ids_by_sequence(sequence):
    """
    Retrieve a list of PDB IDs that match the given sequence.
    """
    get_id_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 0.001,
                "identity_cutoff": 0.9,
                "sequence_type": "protein",
                "value": sequence
            }
        },
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 100
            }
        },
        "return_type": "entry"
    }

    try:
        response = requests.post(get_id_url, json=query, timeout=10)
        if response.status_code == 200:
            results = response.json()
            return [entry['identifier'] for entry in results.get('result_set', [])]
        else:
            print(f"Error fetching PDB IDs: {response.status_code}")
            sys.exit(41)
    except requests.Timeout:
        print(f"Request timed out while fetching PDB IDs.")
        sys.exit(42)
    except requests.RequestException as e:
        print(f"An error occurred while fetching PDB IDs: {e}")
        sys.exit(41)

def get_pdb_info(pdb_id):
    """
    Get information about a PDB entry by its ID.
    """
    get_info_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        response = requests.get(get_info_url, timeout=10)
        if response.status_code == 200:
            result = response.json()
            return result.get('struct', {}).get('title', 'No title available')
        else:
            print(f"Error fetching PDB info for {pdb_id}: {response.status_code}")
            return None
    except requests.Timeout:
        print(f"Request timed out while fetching PDB info.")
        return None
    except requests.RequestException as e:
        print(f"An error occurred while fetching PDB info for {pdb_id}: {e}")
        return None



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
            print("NCBI에 요청한 nucleotide ID가 존재하지 않습니다.")
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
        # Run muscle using subprocess, reading from stdin and writing to stdout
        result = subprocess.run([self.muscle_exe], input=self.combined_memory_file.read(), text=True, capture_output=True)  # 수정된 부분: muscle 실행을 위해 stdin으로 데이터를 전달하고 stdout을 메모리에 저장
        if result.returncode == 0:
            self.aligned_memory_file = StringIO(result.stdout)  # 수정된 부분: 결과를 StringIO 객체에 저장하여 메모리 내에서 처리
        else:
            print(f"Error running MUSCLE: {result.stderr}")
            sys.exit(21)

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
        self.amino_acid_sequence = None
        self.linearDesign = []
        self.protParam = []

    def run_linear_design(self):
        start=self.start-1
        end=self.end

        # Update region
        (idx_start,idx_end) = self.alignment_index[self.gene]
        input_sequence = self.alignment_dict[self.reference_id][idx_start:idx_end]

        # print(input_sequence)
        while True:
            initial_gap_count = input_sequence[:start].count("-")
            gap_count = input_sequence[:start + initial_gap_count].count("-")
            if initial_gap_count == gap_count:
                start += initial_gap_count
                end += initial_gap_count
                break

        while True:
            gap_count_in_range = input_sequence[start:end].count("-")
            updated_end = end + gap_count_in_range
            if updated_end > idx_end:
                end = idx_end
                break
            
            gap_count = input_sequence[start:end + gap_count_in_range].count("-")
            
            if gap_count_in_range == gap_count:
                end += gap_count_in_range
                break
        
        # Get target sequence
        amino_acid_sequence = self.alignment_dict[self.variant_id][idx_start:idx_end][start:end].replace("-", "")
        
        if(amino_acid_sequence == ""):
            print("Error: No amino acid sequence found")
            sys.exit(31)

        # Run LinearDesign
        # Execute the command and capture the result
        try:
            os.chdir(os.path.join(current_dir, "../../LinearDesign"))
        except FileNotFoundError as e:
            print("Error: No LinearDesign directory")
            sys.exit(33)
        command = f"echo {amino_acid_sequence} | ./lineardesign --lambda 3"
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()


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
            self.linearDesign.append(amino_acid_sequence)
            self.linearDesign.append(mRNA_sequence)
            self.linearDesign.append(mRNA_structure)
            self.linearDesign.append(free_energy)
            self.linearDesign.append(cai)
            self.amino_acid_sequence = amino_acid_sequence

        else:
            print(f"Error running LinearDesign: {stderr.decode()}")
            sys.exit(32)
            

    def set_protParam(self):
        # Create a protein analysis object
        protein_analysis = ProteinAnalysis(self.amino_acid_sequence)

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
        self.protParam.append(molecular_weight)
        self.protParam.append(amino_acid_count)
        self.protParam.append(amino_acid_percent)
        self.protParam.append(isoelectric_point)
        self.protParam.append(instability_index)
        self.protParam.append(secondary_structure_fraction)
        self.protParam.append(gravy)
        self.protParam.append(aromaticity)

    def get_linearDesign(self):
        amino_acid_sequence, mRNA_sequence, mRNA_structure, free_energy, cai= self.linearDesign
        linearDesign_dict = {
            "amino_acid_sequence": amino_acid_sequence,
            "mRNA_sequence": mRNA_sequence,
            "mRNA_structure": mRNA_structure,
            "free_energy": free_energy,
            "cai": cai
        }
        return linearDesign_dict

    def get_protParam(self):
        molecular_weight, amino_acid_count, amino_acid_percent, isoelectric_point, instability_index, secondary_structure_fraction, gravy, aromaticity = self.protParam
        protParam_dict = {
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
    current_dir=sys.argv[1]
    option = int(sys.argv[2])

    # metadata
    if option == 1:
        reference_id = sys.argv[3]
        os.makedirs(current_dir+"/data", exist_ok=True)
        metadata = get_metadata(reference_id)
        save_json(metadata, "metadata.json")  # JSON 파일로 저장

    # alignment
    elif option == 2:
        variant_sequences = {}
        # read fasta file
        fasta_content = sys.argv[3]
        fasta_io = StringIO(fasta_content)
        for record in SeqIO.parse(fasta_io, "fasta"):
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

        save_json(alignment_data, "alignment.json")  # JSON 파일로 저장

    # linearDesign, protparam data
    elif option == 3:
        # get metadata
        metadata = get_json("metadata.json")
        reference_id = metadata.get("Sequence ID", None)
        
        # get alignment data
        alignment_data = get_json("alignment.json")
        alignment_index = alignment_data.get("alignment_index", None)
        alignment_dict = alignment_data.get("aligned_sequences", None)

        # set gene, variant_id, start, end
        gene=sys.argv[3]
        variant_id=sys.argv[4]
        start = int(sys.argv[5])
        end = int(sys.argv[6])

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
        save_json(linearDesign_data, "linearDesign.json")  # JSON 파일로 저장


    elif option == 4:
        # get metadata
        metadata = get_json("metadata.json")
        reference_id = metadata.get("Sequence ID", None)
        
        # get alignment data
        alignment_data = get_json("alignment.json")
        alignment_index = alignment_data.get("alignment_index", None)
        alignment_dict = alignment_data.get("aligned_sequences", None)

        # set gene, variant_id, start, end
        gene=sys.argv[3]
        
        # PDB search
        sequence = alignment_dict[reference_id][alignment_index[gene][0]:alignment_index[gene][1]].replace("-", "")

        pdb_ids = get_pdb_ids_by_sequence(sequence)
        pdb_ids=pdb_ids[:10]
        pdb_dict = {pdb_id: "" for pdb_id in pdb_ids}

        for pdb_id in pdb_ids:
            title = get_pdb_info(pdb_id)
            if title:
                pdb_dict[pdb_id] = title
            else:
                break

        save_json(pdb_dict, "pdb.json")