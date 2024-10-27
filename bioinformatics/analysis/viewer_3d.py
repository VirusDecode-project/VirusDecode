import json
import requests
import sys

class ThreeDViewer:
    def __init__(self, sequence):
        self.sequence = sequence
        self.pdb_ids = []
        self.pdb_dict = {}

    def get_pdb_ids_by_sequence(self):
        """Retrieve a list of PDB IDs that match the given sequence."""
        get_id_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        query = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": 0.001,
                    "identity_cutoff": 0.9,
                    "sequence_type": "protein",
                    "value": self.sequence
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
                self.pdb_ids = [entry['identifier'] for entry in results.get('result_set', [])][:10]
            else:
                sys.stderr.write(f"Error fetching PDB IDs: {response.status_code}\n")
                sys.exit(41)
        except requests.Timeout:
            sys.stderr.write("Request timed out while fetching PDB IDs.\n")
            sys.exit(42)
        except requests.RequestException as e:
            sys.stderr.write(f"An error occurred while fetching PDB IDs: {str(e)}\n")
            sys.exit(41)

    def get_pdb_info(self, pdb_id):
        """Get information about a PDB entry by its ID."""
        get_info_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        
        try:
            response = requests.get(get_info_url, timeout=10)
            if response.status_code == 200:
                result = response.json()
                return result.get('struct', {}).get('title', 'No title available')
            else:
                return None
        except requests.Timeout:
            return None
        except requests.RequestException:
            return None

    def fetch_pdb_info(self):
        """Fetch titles for all retrieved PDB IDs and store them in a dictionary."""
        self.pdb_dict = {pdb_id: "" for pdb_id in self.pdb_ids}
        for pdb_id in self.pdb_ids:
            title = self.get_pdb_info(pdb_id)
            if title:
                self.pdb_dict[pdb_id] = title
            else:
                break

    def print_pdb_info(self):
        """Print the PDB information as a JSON-formatted string."""
        print(json.dumps(self.pdb_dict, indent=4))

    def run(self):
        """Run the entire fetch and print process."""
        self.get_pdb_ids_by_sequence()
        self.fetch_pdb_info()
        self.print_pdb_info()