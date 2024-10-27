import unittest
from unittest.mock import patch, MagicMock, mock_open
from io import StringIO
import json
from Bio.Seq import Seq
from virusdecode import get_metadata, get_pdb_ids_by_sequence, get_pdb_info, SequenceAlignment, SequenceAnalysis, is_valid_sequence

class TestFunctions(unittest.TestCase):

    @patch("virusdecode.Entrez.efetch")
    @patch("virusdecode.SeqIO.read")
    def test_get_metadata(self, mock_read, mock_efetch):
        # Mock Entrez response
        mock_handle = MagicMock()
        mock_efetch.return_value = mock_handle
        mock_reference = MagicMock()
        mock_reference.id = "TestID"
        mock_reference.name = "TestName"
        mock_reference.description = "TestDescription"
        mock_reference.__len__.return_value = 100
        mock_read.return_value = mock_reference

        result = get_metadata("dummy_id")
        expected = {
            "Sequence ID": "TestID",
            "Name": "TestName",
            "Description": "TestDescription",
            "Length": 100
        }
        self.assertEqual(result, expected)

    @patch("virusdecode.requests.post")
    def test_get_pdb_ids_by_sequence(self, mock_post):
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "result_set": [{"identifier": "1XYZ"}, {"identifier": "2ABC"}]
        }
        mock_post.return_value = mock_response

        result = get_pdb_ids_by_sequence("dummy_sequence")
        expected = ["1XYZ", "2ABC"]
        self.assertEqual(result, expected)

    @patch("virusdecode.requests.get")
    def test_get_pdb_info(self, mock_get):
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"struct": {"title": "Test Title"}}
        mock_get.return_value = mock_response

        result = get_pdb_info("dummy_pdb_id")
        self.assertEqual(result, "Test Title")

    def test_is_valid_sequence(self):
        self.assertTrue(is_valid_sequence("ATCGATCG"))
        self.assertFalse(is_valid_sequence("AXTG"))

class TestSequenceAlignment(unittest.TestCase):

    @patch("virusdecode.Entrez.efetch")
    @patch("virusdecode.SeqIO.read")
    def test_sequence_alignment_initialization(self, mock_read, mock_efetch):
        mock_handle = MagicMock()
        mock_efetch.return_value = mock_handle
        mock_reference = MagicMock()
        mock_reference.id = "TestID"
        mock_read.return_value = mock_reference

        alignment = SequenceAlignment({"var1": Seq("ATCG")}, "NC_045512")
        self.assertEqual(alignment.reference_id, "TestID")
        self.assertIn("var1", alignment.variant_sequences)

    @patch("virusdecode.subprocess.run")
    def test_run_muscle_dna(self, mock_run):
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = ">TestID\nATCG\n>var1\nTAGC\n"
        
        alignment = SequenceAlignment({"var1": Seq("ATCG")}, "NC_045512")
        alignment.combined_memory_file = StringIO(">TestID\nATCG\n>var1\nTAGC\n")
        alignment.run_muscle_dna()
        result = alignment.aligned_memory_file.getvalue()
        
        self.assertIn("ATCG", result)
        self.assertIn("TAGC", result)

class TestSequenceAnalysis(unittest.TestCase):

    @patch("virusdecode.subprocess.Popen")
    def test_run_linear_design(self, mock_popen):
        # `subprocess.Popen`의 `communicate` 메서드가 반환할 값 설정
        mock_process = MagicMock()
        mock_process.communicate.return_value = (
            b"mRNA sequence:  AUGUUCGUGUUCCUAGUGCUGCUGCCACUAGUGAGCUCG\nmRNA structure: .....((.(((((((((((....).)))))).)))).))\nmRNA folding free energy: -13.40 kcal/mol; mRNA CAI: 0.674" , 
            b""
        )
        mock_process.returncode = 0
        mock_popen.return_value = mock_process

        analysis = SequenceAnalysis("MFVFLVLLPLVSS")
        analysis.run_linear_design()  # 예상된 결과로 실행
        result = analysis.get_linearDesign()

        expected = {
            "amino_acid_sequence": "MFVFLVLLPLVSS",
            "mRNA_sequence": "AUGUUCGUGUUCCUAGUGCUGCUGCCACUAGUGAGCUCG",
            "mRNA_structure": ".....((.(((((((((((....).)))))).)))).))",
            "free_energy": "-13.40 kcal/mol",
            "cai": "0.674"
        }
        self.assertEqual(result, expected)

    @patch("virusdecode.ProteinAnalysis")
    def test_set_protParam(self, mock_protein_analysis):
        mock_protein = MagicMock()
        mock_protein.molecular_weight.return_value = 2820.4148000000005
        mock_protein.count_amino_acids.return_value = {'A': 0, 'C': 1, 'D': 0, 'E': 0, 'F': 2, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 6, 'M': 1, 'N': 1, 'P': 2, 'Q': 2, 'R': 1, 'S': 2, 'T': 3, 'V': 4, 'W': 0, 'Y': 0}
        mock_protein.get_amino_acids_percent.return_value = {'A': 0.0, 'C': 0.04, 'D': 0.0, 'E': 0.0, 'F': 0.08, 'G': 0.0, 'H': 0.0, 'I': 0.0, 'K': 0.0, 'L': 0.24, 'M': 0.04, 'N': 0.04, 'P': 0.08, 'Q': 0.08, 'R': 0.04, 'S': 0.08, 'T': 0.12, 'V': 0.16, 'W': 0.0, 'Y': 0.0}
        mock_protein.isoelectric_point.return_value = 7.999833488464356
        mock_protein.instability_index.return_value = 31.368000000000002
        mock_protein.secondary_structure_fraction.return_value = (0.27999999999999997, 0.2, 0.6)
        mock_protein.gravy.return_value = 1.108
        mock_protein.aromaticity.return_value = 0.08
        mock_protein_analysis.return_value = mock_protein

        analysis = SequenceAnalysis("MESLVPGFNEKTHVQLSLPVLQVRD")
        analysis.set_protParam()
        result = analysis.get_protParam()
        
        expected = {
            'molecular_weight': 2820.4148000000005, 
            'amino_acid_count': {'A': 0, 'C': 1, 'D': 0, 'E': 0, 'F': 2, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 6, 'M': 1, 'N': 1, 'P': 2, 'Q': 2, 'R': 1, 'S': 2, 'T': 3, 'V': 4, 'W': 0, 'Y': 0}, 
            'amino_acid_percent': {'A': 0.0, 'C': 0.04, 'D': 0.0, 'E': 0.0, 'F': 0.08, 'G': 0.0, 'H': 0.0, 'I': 0.0, 'K': 0.0, 'L': 0.24, 'M': 0.04, 'N': 0.04, 'P': 0.08, 'Q': 0.08, 'R': 0.04, 'S': 0.08, 'T': 0.12, 'V': 0.16, 'W': 0.0, 'Y': 0.0}, 
            'isoelectric_point': 7.999833488464356, 
            'instability_index': 31.368000000000002, 
            'secondary_structure_fraction': (0.27999999999999997, 0.2, 0.6), 
            'gravy': 1.108, 
            'aromaticity': 0.08
        }
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
