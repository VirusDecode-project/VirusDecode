import unittest
from unittest.mock import patch, MagicMock, mock_open
from io import StringIO
import sys
import json
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import subprocess
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../analysis')))
from mRNA_design import MRNADesign

class TestMRNADesignSuccessCases(unittest.TestCase):
    
    @patch("mRNA_design.subprocess.Popen")
    def test_run_linear_design_success(self, mock_popen):
        mock_process = MagicMock()
        mock_process.communicate.return_value = (
            b"mRNA sequence: AUGUUCGUGUUCCUAGUGCUGCUGCCACUAGUGAGCUCG\nmRNA structure: .....((.(((((((((((....).)))))).)))).))\nmRNA folding free energy: -13.40 kcal/mol; mRNA CAI: 0.674",
            b""
        )
        mock_process.returncode = 0
        mock_popen.return_value = mock_process

        analysis = MRNADesign("MFVFLVLLPLVSS")
        analysis.run_linear_design()
        result = analysis.get_linearDesign()

        expected = {
            "amino_acid_sequence": "MFVFLVLLPLVSS",
            "mRNA_sequence": "AUGUUCGUGUUCCUAGUGCUGCUGCCACUAGUGAGCUCG",
            "mRNA_structure": ".....((.(((((((((((....).)))))).)))).))",
            "free_energy": "-13.40 kcal/mol",
            "cai": "0.674"
        }
        self.assertEqual(result, expected)

    @patch("mRNA_design.ProteinAnalysis")
    def test_set_protParam_success(self, mock_protein_analysis):
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

        analysis = MRNADesign("MESLVPGFNEKTHVQLSLPVLQVRD")
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

class TestMRNADesignFailureCases(unittest.TestCase):
    
    @patch("mRNA_design.subprocess.Popen")
    @patch("mRNA_design.sys.exit")
    @patch("mRNA_design.sys.stderr", new_callable=StringIO)
    def test_run_linear_design_failure(self, mock_stderr, mock_exit, mock_popen):
        # Simulate a failure in the LinearDesign command
        process_mock = MagicMock()
        mock_popen.return_value = process_mock
        process_mock.communicate.return_value = ("", "Error occurred")
        process_mock.returncode = 1

        mrna_design = MRNADesign("MESLVPGFNEKTHVQLSLPVLQVRD")
        mrna_design.run_linear_design()

        self.assertIn("Command failed with return code", mock_stderr.getvalue())
        mock_exit.assert_called_once_with(32)

    @patch("mRNA_design.ProteinAnalysis")
    def test_set_protParam_invalid_sequence(self, mock_protein_analysis):
        # Test an invalid amino acid sequence
        mock_protein_analysis.side_effect = ValueError("Invalid sequence")

        mrna_design = MRNADesign("INVALID")
        with self.assertRaises(ValueError):
            mrna_design.set_protParam()

    @patch("mRNA_design.os.chdir", side_effect=FileNotFoundError("No such directory"))
    @patch("mRNA_design.sys.exit")
    @patch("mRNA_design.sys.stderr", new_callable=StringIO)
    def test_run_linear_design_directory_not_found(self, mock_stderr, mock_exit, mock_chdir):
        mrna_design = MRNADesign("MESLVPGFNEKTHVQLSLPVLQVRD")
        mrna_design.run_linear_design()

        self.assertIn("Not exist LinearDesign directory.", mock_stderr.getvalue())
        mock_exit.assert_called_once_with(33)

if __name__ == "__main__":
    unittest.main()