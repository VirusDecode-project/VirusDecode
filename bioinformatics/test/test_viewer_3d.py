import unittest
from unittest.mock import patch, MagicMock
import json
from io import StringIO
import sys
import os
import requests
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../analysis')))
from viewer_3d import ThreeDViewer

class TestThreeDViewerSuccessCases(unittest.TestCase):

    @patch("viewer_3d.requests.post")
    def test_get_pdb_ids_by_sequence_success(self, mock_post):
        # Mock the response of the requests.post for get_pdb_ids_by_sequence
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "result_set": [
                {"identifier": "1ABC"},
                {"identifier": "2DEF"},
                {"identifier": "3GHI"}
            ]
        }
        mock_post.return_value = mock_response

        viewer = ThreeDViewer("MESLVPGFNEKTHVQLSLPVLQVRD")
        viewer.get_pdb_ids_by_sequence()

        self.assertEqual(viewer.pdb_ids, ["1ABC", "2DEF", "3GHI"])

    @patch("viewer_3d.requests.get")
    def test_get_pdb_info_success(self, mock_get):
        # Mock the response of requests.get for get_pdb_info
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "struct": {
                "title": "Test Protein Structure"
            }
        }
        mock_get.return_value = mock_response

        viewer = ThreeDViewer("MESLVPGFNEKTHVQLSLPVLQVRD")
        title = viewer.get_pdb_info("1ABC")

        self.assertEqual(title, "Test Protein Structure")

    @patch("viewer_3d.ThreeDViewer.get_pdb_info", return_value="Sample Protein")
    def test_fetch_pdb_info_success(self, mock_get_pdb_info):
        viewer = ThreeDViewer("MESLVPGFNEKTHVQLSLPVLQVRD")
        viewer.pdb_ids = ["1ABC", "2DEF"]
        viewer.fetch_pdb_info()

        expected_pdb_dict = {
            "1ABC": "Sample Protein",
            "2DEF": "Sample Protein"
        }
        self.assertEqual(viewer.pdb_dict, expected_pdb_dict)

    @patch("viewer_3d.ThreeDViewer.print_pdb_info")
    @patch("viewer_3d.ThreeDViewer.fetch_pdb_info")
    @patch("viewer_3d.ThreeDViewer.get_pdb_ids_by_sequence")
    def test_run_success(self, mock_get_pdb_ids, mock_fetch_pdb_info, mock_print_pdb_info):
        viewer = ThreeDViewer("MESLVPGFNEKTHVQLSLPVLQVRD")
        viewer.run()

        mock_get_pdb_ids.assert_called_once()
        mock_fetch_pdb_info.assert_called_once()
        mock_print_pdb_info.assert_called_once()

class TestThreeDViewerFailureCases(unittest.TestCase):

    @patch("viewer_3d.requests.post", side_effect=requests.Timeout)
    @patch("viewer_3d.sys.stderr", new_callable=StringIO)
    @patch("viewer_3d.sys.exit")
    def test_get_pdb_ids_by_sequence_timeout(self, mock_exit, mock_stderr, mock_post):
        viewer = ThreeDViewer("MESLVPGFNEKTHVQLSLPVLQVRD")
        viewer.get_pdb_ids_by_sequence()

        self.assertIn("Request timed out while fetching PDB IDs", mock_stderr.getvalue())
        mock_exit.assert_called_once_with(42)

    @patch("viewer_3d.requests.get", side_effect=requests.Timeout)
    def test_get_pdb_info_timeout(self, mock_get):
        viewer = ThreeDViewer("MESLVPGFNEKTHVQLSLPVLQVRD")
        title = viewer.get_pdb_info("1ABC")

        self.assertIsNone(title)

    @patch("viewer_3d.requests.get", side_effect=requests.RequestException("Connection error"))
    def test_get_pdb_info_request_exception(self, mock_get):
        viewer = ThreeDViewer("MESLVPGFNEKTHVQLSLPVLQVRD")
        title = viewer.get_pdb_info("1ABC")

        self.assertIsNone(title)

if __name__ == "__main__":
    unittest.main()
