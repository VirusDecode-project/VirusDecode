import requests

# RCSB PDB API endpoint
url = "https://search.rcsb.org/rcsbsearch/v2/query"

NCBI_id = "NC_045512.2" # CoV2 RNA

query = {
    "query": {
        "type": "terminal",
        "service": "full_text",
        "parameters": {
            "value": NCBI_id
        }
    },
    "return_type": "entry"
}

# request API
response = requests.post(url, json=query)

# get PDB IDs (list)
if response.status_code == 200:
    results = response.json()
    pdb_ids = [entry['identifier'] for entry in results['result_set']]
    print("PDB IDs:", pdb_ids[:10]) # max pdb id number : 10
else:
    print(f"Error: {response.status_code}")

