import json
import sys
from Bio import SeqIO, Entrez
Entrez.email = "your_email@example.com"

class Metadata:
    def __init__(self, reference_id):
        self.reference_id = reference_id
        self.metadata = {}
    def set_metadata(self):
        try:
            handle = Entrez.efetch(db="nucleotide", id=self.reference_id, rettype="gb", retmode="text")
            reference_sequence = SeqIO.read(handle, "genbank")
            handle.close()
            metadata={
                "Sequence ID": reference_sequence.id,
                "Name": reference_sequence.name,
                "Description": reference_sequence.description,
                "Length": len(reference_sequence)
            }
            self.metadata = metadata
        except Exception as e:
            sys.stderr.write(f"There was an error fetching the reference sequence from NCBI.\n{str(e)}")
            sys.exit(11)

    def print_metadata(self):
        print(json.dumps(self.metadata, indent=4))

    def run(self):
        self.set_metadata()
        self.print_metadata()