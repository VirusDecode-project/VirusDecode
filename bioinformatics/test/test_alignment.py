import unittest
from unittest.mock import patch, MagicMock, mock_open
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../analysis')))
from alignment import Alignment

class TestAlignmentSuccessCases(unittest.TestCase):
    
    @patch("alignment.Entrez.efetch")
    @patch("alignment.SeqIO.read")
    def test_validate_and_fetch_reference_sequence_success(self, mock_seqio_read, mock_efetch):
        """성공적으로 레퍼런스 시퀀스를 가져와서 저장하는지 확인"""
        mock_handle = MagicMock()
        mock_efetch.return_value = mock_handle
        mock_reference_sequence = SeqRecord(Seq("ATCG"))
        mock_reference_sequence.id = "TestID"
        mock_reference_sequence.features = []
        mock_seqio_read.return_value = mock_reference_sequence

        alignment = Alignment("dummy_id", "fasta_content")
        alignment.fetch_reference_sequence()

        self.assertEqual(alignment.reference_sequence.id, "TestID")

    @patch("alignment.subprocess.run")
    @patch("alignment.Entrez.efetch")
    @patch("alignment.SeqIO.read")
    def test_run_muscle_dna_success(self, mock_seqio_read, mock_efetch, mock_subprocess_run):
        """MUSCLE 실행 후 결과를 메모리에 저장하는지 확인"""
        mock_handle = MagicMock()
        mock_efetch.return_value = mock_handle
        mock_reference_sequence = SeqRecord(Seq("ATCG"))
        mock_reference_sequence.id = "TestID"
        mock_reference_sequence.features = []
        mock_seqio_read.return_value = mock_reference_sequence
        mock_subprocess_run.return_value = MagicMock(stdout=">TestID\nATCG\n>Variant1\nTAGCAA")

        alignment = Alignment("dummy_id", ">Variant1\nTAGCAA")
        alignment.run()

        # 결과가 예상대로 저장되었는지 확인
        self.assertEqual(alignment.aligned_memory_file.getvalue(), ">TestID\nATCG\n>Variant1\nTAGCAA")

class TestAlignmentFailureCases(unittest.TestCase):

    @patch("alignment.Entrez.efetch", side_effect=Exception("Fetch failed"))
    @patch("alignment.sys.exit")
    @patch("alignment.sys.stderr", new_callable=StringIO)
    def test_fetch_reference_sequence_failure(self, mock_stderr, mock_exit, mock_efetch):
        """fetch_reference_sequence에서 예외 발생 시 처리 확인"""
        alignment = Alignment("dummy_id", "fasta_content")
        alignment.fetch_reference_sequence()

        # 에러 메시지 확인
        self.assertIn("There was an error fetching the reference sequence from NCBI.", mock_stderr.getvalue())
        mock_exit.assert_called_once_with(11)

    @patch("alignment.SeqIO.parse", side_effect=ValueError("Parse failed"))
    @patch("alignment.sys.exit")
    @patch("alignment.sys.stderr", new_callable=StringIO)
    def test_validate_sequences_parsing_error(self, mock_stderr, mock_exit, mock_seqio_parse):
        """validate_sequences에서 파싱 오류 발생 시 처리 확인"""
        alignment = Alignment("dummy_id", "invalid_fasta_content")
        alignment.validate_sequences()

        # stderr에 파싱 오류 메시지 확인
        self.assertIn("fasta format error in SeqIO.parse.", mock_stderr.getvalue())
        mock_exit.assert_called_once_with(22)

    @patch("alignment.subprocess.run", side_effect=subprocess.CalledProcessError(1, 'muscle'))
    def test_run_muscle_dna_failure(self, mock_run):
        alignment = Alignment({"var1": Seq("ATCG")}, "NC_045512")
        alignment.combined_memory_file = StringIO(">TestID\nATCG\n>var1\nTAGC\n")

        with self.assertRaises(SystemExit) as cm:
            with patch('sys.stderr', new=StringIO()) as mock_stderr:
                alignment.run_muscle_dna()
                self.assertIn("Error running MUSCLE", mock_stderr.getvalue())
        self.assertEqual(cm.exception.code, 21)

if __name__ == "__main__":
    unittest.main()
