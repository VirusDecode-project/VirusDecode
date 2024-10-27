import unittest
from unittest.mock import patch, MagicMock
import json
import sys
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../analysis')))
from metadata import Metadata

class TestMetadataSuccessCases(unittest.TestCase):
    
    @patch("metadata.Entrez.efetch")
    @patch("metadata.SeqIO.read")
    def test_set_metadata_success(self, mock_seqio_read, mock_efetch):
        # Mocking Entrez handle and SeqRecord response
        mock_handle = MagicMock()
        mock_efetch.return_value = mock_handle
        mock_reference_sequence = SeqRecord(Seq("ATCG"))
        mock_reference_sequence.id = "TestID"
        mock_reference_sequence.name = "TestName"
        mock_reference_sequence.description = "TestDescription"
        mock_seqio_read.return_value = mock_reference_sequence

        metadata = Metadata("dummy_id")
        metadata.set_metadata()

        expected_metadata = {
            "Sequence ID": "TestID",
            "Name": "TestName",
            "Description": "TestDescription",
            "Length": 4
        }
        self.assertEqual(metadata.metadata, expected_metadata)

    @patch("metadata.sys.stderr", new_callable=StringIO)
    @patch("metadata.sys.exit")
    @patch("metadata.Entrez.efetch", side_effect=Exception("Fetch failed"))
    def test_set_metadata_http_error(self, mock_efetch, mock_exit, mock_stderr):
        metadata = Metadata("dummy_id")
        
        metadata.set_metadata()
        
        # Check that error message was written to stderr and sys.exit was called
        self.assertIn("There was an error fetching the reference sequence from NCBI.", mock_stderr.getvalue())
        mock_exit.assert_called_once_with(11)

    @patch("metadata.Metadata.set_metadata")
    @patch("metadata.sys.stdout", new_callable=StringIO)
    def test_print_metadata(self, mock_stdout, mock_set_metadata):
        # Prepare mock metadata
        mock_set_metadata.return_value = None
        metadata = Metadata("dummy_id")
        metadata.metadata = {
            "Sequence ID": "TestID",
            "Name": "TestName",
            "Description": "TestDescription",
            "Length": 4
        }
        
        metadata.print_metadata()
        
        expected_output = json.dumps(metadata.metadata, indent=4) + "\n"
        self.assertEqual(mock_stdout.getvalue(), expected_output)

    @patch("metadata.Metadata.set_metadata")
    @patch("metadata.Metadata.print_metadata")
    def test_run(self, mock_print_metadata, mock_set_metadata):
        metadata = Metadata("dummy_id")
        metadata.run()

        # Ensure set_metadata and print_metadata were called in run
        mock_set_metadata.assert_called_once()
        mock_print_metadata.assert_called_once()

class TestMetadataFailureCases(unittest.TestCase):
    
    @patch("metadata.Entrez.efetch", side_effect=Exception("Fetch failed"))
    @patch("metadata.sys.exit")
    @patch("metadata.sys.stderr", new_callable=StringIO)
    def test_set_metadata_invalid_reference_id(self, mock_stderr, mock_exit, mock_efetch):
        """잘못된 reference_id로 Exception가 발생하는 경우 테스트"""
        metadata = Metadata("invalid_id")
        metadata.set_metadata()
        
        self.assertIn("There was an error fetching the reference sequence from NCBI.", mock_stderr.getvalue())
        mock_exit.assert_called_once_with(11)

    @patch("metadata.Entrez.efetch")
    @patch("metadata.SeqIO.read", side_effect=Exception("Fetch failed"))
    @patch("metadata.sys.exit")
    @patch("metadata.sys.stderr", new_callable=StringIO)
    def test_set_metadata_parsing_error(self, mock_stderr, mock_exit, mock_seqio_read, mock_efetch):
        """SeqIO가 시퀀스를 파싱하지 못해 ValueError가 발생하는 경우 테스트"""
        metadata = Metadata("dummy_id")
        metadata.set_metadata()

        self.assertIn("There was an error fetching the reference sequence from NCBI.", mock_stderr.getvalue())
        mock_exit.assert_called_once_with(11)

    @patch("metadata.Entrez.efetch")
    @patch("metadata.SeqIO.read")
    def test_set_metadata_incomplete_data(self, mock_seqio_read, mock_efetch):
        """레퍼런스 데이터가 불완전한 경우 테스트 (ID 없음)"""
        mock_handle = MagicMock()
        mock_efetch.return_value = mock_handle
        mock_reference_sequence = SeqRecord(Seq("ATCG"))
        mock_reference_sequence.id = None  # ID가 없는 경우
        mock_reference_sequence.name = "TestName"
        mock_reference_sequence.description = "TestDescription"
        mock_seqio_read.return_value = mock_reference_sequence

        metadata = Metadata("dummy_id")
        metadata.set_metadata()

        self.assertEqual(metadata.metadata.get("Sequence ID"), None)
        self.assertEqual(metadata.metadata["Name"], "TestName")
        self.assertEqual(metadata.metadata["Description"], "TestDescription")

if __name__ == "__main__":
    unittest.main()

