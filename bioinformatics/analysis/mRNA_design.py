import os
import sys
import json
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import subprocess
    
class MRNADesign:
    def __init__(self, amino_acid_sequence):
        self.amino_acid_sequence = amino_acid_sequence
        self.linearDesign = []
        self.protParam = []
        self.linearDesign_data = {}

    def run_linear_design(self):
        # Run LinearDesign
        try:
            os.chdir(os.path.join(os.getcwd(), "../../LinearDesign"))
            command = f"echo {self.amino_acid_sequence} | ./lineardesign --lambda 3"
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            # Check if the command was executed successfully
            if process.returncode != 0:
                sys.stderr.write(f"Command failed with return code {process.returncode}.\n")
                sys.exit(32)
            else:
                if not stdout.decode().strip():
                    sys.stderr.write("Error in running time by LinearDesign. It could be timeout.\n")
                    sys.exit(31)

                # Parse the output as needed
                output_lines = [line for line in stdout.decode().splitlines() if not line.startswith("j=")]
                if len(output_lines) >= 3:
                    mRNA_sequence = output_lines[0].replace('mRNA sequence:', '').strip()
                    mRNA_structure = output_lines[1].replace('mRNA structure:', '').strip()
                    parts = output_lines[2].split(';')
                    free_energy = parts[0].replace('mRNA folding free energy:', '').strip()
                    cai = parts[1].replace('mRNA CAI:', '').strip()

                    # Set the linear design data
                    self.linearDesign.append(self.amino_acid_sequence)
                    self.linearDesign.append(mRNA_sequence)
                    self.linearDesign.append(mRNA_structure)
                    self.linearDesign.append(free_energy)
                    self.linearDesign.append(cai)
                else:
                    sys.stderr.write("Error in running time by LinearDesign. It could be timeout.\n")
                    sys.exit(31)
        except FileNotFoundError as e:
            sys.stderr.write(f"Not exist LinearDesign directory.\n{str(e)}")
            sys.exit(33)

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


    def set_linearDesign_data(self):
        # get linearDesign and protParam data
        linearDesign_dict = self.get_linearDesign()
        protParam_dict = self.get_protParam()
        linearDesign_data = {
            "linearDesign": linearDesign_dict,
            "protParam": protParam_dict
        }

        self.linearDesign_data = linearDesign_data

    def print_linearDesign_data(self):
        print(json.dumps(self.linearDesign_data, indent=4))

    def run(self):
        self.run_linear_design()
        self.set_protParam()
        self.set_linearDesign_data()
        self.print_linearDesign_data()