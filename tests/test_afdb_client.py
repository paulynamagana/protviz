import logging
import unittest
from unittest.mock import Mock, call, patch

import requests

# Assuming the AFDBClient is in this path for testing
from protviz.data_retrieval.afdb_client import AFDBClient

# --- Sample data for mocking (unchanged) ---
SAMPLE_UNIPROT_ID = "P12345"
SAMPLE_API_RESPONSE_ENTRY = [
    {
        "uniprotAccession": SAMPLE_UNIPROT_ID,
        "uniprotId": "TEST_ID",
        "gene": "TEST_GENE",
        "latestVersion": 3,
        "cifUrl": f"https://alphafold.ebi.ac.uk/files/AF-{SAMPLE_UNIPROT_ID}-F1-model_v3.cif",
        "pdbUrl": f"https://alphafold.ebi.ac.uk/files/AF-{SAMPLE_UNIPROT_ID}-F1-model_v3.pdb",
        "paeImageUrl": f"https://alphafold.ebi.ac.uk/files/AF-{SAMPLE_UNIPROT_ID}-F1-pae_v3.png",
        "amAnnotationsUrl": f"https://alphafold.ebi.ac.uk/files/AF-{SAMPLE_UNIPROT_ID}-F1-alphafold_predictions_v3.csv",
        "taxId": 9606,  # Human
    }
]
SAMPLE_CIF_CONTENT = """
#
data_AF-P12345-F1-modelv3
#
loop_
_ma_qa_metric_local.label_seq_id
_ma_qa_metric_local.metric_value
1 90.5
2 85.3
3 95.0
#
"""
SAMPLE_CIF_CONTENT_MISMATCHED_LENGTHS = """
#
data_AF-P12345-F1-modelv3
#
loop_
_ma_qa_metric_local.label_seq_id
_ma_qa_metric_local.metric_value
1 90.5
2
3 95.0
#
"""
SAMPLE_CIF_CONTENT_BAD_VALUES = """
#
data_AF-P12345-F1-modelv3
#
loop_
_ma_qa_metric_local.label_seq_id
_ma_qa_metric_local.metric_value
A 90.5
2 B
3 95.0
#
"""
SAMPLE_ALPHAMISSENSE_CSV_CONTENT = """protein_variant,am_pathogenicity,am_class
M1A,0.123,likely_benign
L2P,0.890,likely_pathogenic
P3X,0.500,ambiguous
"""
SAMPLE_ALPHAMISSENSE_CSV_CONTENT_BAD_HEADER = """variant,score,classification
M1A,0.123,likely_benign
"""
SAMPLE_ALPHAMISSENSE_CSV_CONTENT_SHORT_ROW = """protein_variant,am_pathogenicity,am_class
M1A,0.123
"""
SAMPLE_ALPHAMISSENSE_CSV_CONTENT_INVALID_VARIANT = """protein_variant,am_pathogenicity,am_class
INVALID,0.123,likely_benign
M123BADC,0.5,likely_benign
"""


class TestAFDBClient(unittest.TestCase):
    """
    Updated test suite for the refactored AFDBClient.
    NOTE: This test suite assumes that the method names in the source script
    have been corrected for consistency (e.g., `init` -> `__init__`).
    """

    def setUp(self):
        """Set up a test client with patched dependencies for network and filesystem."""
        # Define module path for clarity
        client_module_path = "protviz.data_retrieval.afdb_client"

        # Patch dependencies where they are *used* to prevent side effects
        self.platformdirs_patcher = patch(
            f"{client_module_path}.platformdirs.user_cache_dir",
            return_value="/fake/cache/dir",
        )
        # *** ADDED PATCH FOR os.makedirs TO FIX THE OSError ***
        self.makedirs_patcher = patch(f"{client_module_path}.os.makedirs")
        self.requests_cache_patcher = patch(
            f"{client_module_path}.requests_cache.CachedSession"
        )

        # Start the patchers
        self.mock_user_cache_dir = self.platformdirs_patcher.start()
        self.mock_makedirs = self.makedirs_patcher.start()
        self.mock_cached_session_cls = self.requests_cache_patcher.start()

        # The AFDBClient can now be safely instantiated without side-effects.
        self.client = AFDBClient(timeout=5)
        # Get the instance of the mock session for configuration in tests
        self.mock_session = self.client.session

        # Suppress logging for clean test output
        logging.disable(logging.CRITICAL)

    def tearDown(self):
        """Stop all patchers and re-enable logging."""
        self.platformdirs_patcher.stop()
        self.makedirs_patcher.stop()
        self.requests_cache_patcher.stop()
        logging.disable(logging.NOTSET)

    # ===============================================================
    # The rest of the test methods remain the same as the previous update
    # ===============================================================

    def test_make_api_request_success(self):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"data": "success"}
        mock_response.content = b'{"data": "success"}'
        self.mock_session.get.return_value = mock_response

        # Assumes method is named `_make_api_request`
        response = self.client._make_api_request("/test_endpoint")

        self.assertEqual(response, {"data": "success"})
        self.mock_session.get.assert_called_once_with(
            f"{self.client.AFDB_API_BASE_URL}/test_endpoint",
            timeout=self.client.timeout,
            headers={"Accept": "application/json"},
        )
        mock_response.raise_for_status.assert_called_once()

    def test_make_api_request_empty_response(self):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b""
        self.mock_session.get.return_value = mock_response
        response = self.client._make_api_request("/empty_endpoint")
        self.assertIsNone(response)
        mock_response.raise_for_status.assert_called_once()

    def test_make_api_request_http_error(self):
        mock_response = Mock()
        mock_response.status_code = 404
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
            response=mock_response
        )
        self.mock_session.get.return_value = mock_response
        with self.assertRaises(requests.exceptions.HTTPError):
            self.client._make_api_request("/not_found")

    def test_make_api_request_json_decode_error(self):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b"not a json"
        mock_response.text = "not a json"
        mock_response.json.side_effect = requests.exceptions.JSONDecodeError(
            "msg", "doc", 0
        )
        self.mock_session.get.return_value = mock_response
        with self.assertRaises(ValueError) as cm:
            self.client._make_api_request("/invalid_json")
        self.assertIn("Invalid JSON response", str(cm.exception))

    def test_make_api_request_request_exception(self):
        self.mock_session.get.side_effect = requests.exceptions.RequestException(
            "Network error"
        )
        with self.assertRaises(requests.exceptions.RequestException):
            self.client._make_api_request("/network_error")

    def test_fetch_file_content_success(self):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b"file content"
        self.mock_session.get.return_value = mock_response

        # Assumes method is named `_fetch_file_content`
        content = self.client._fetch_file_content("http://example.com/file.txt")

        self.assertEqual(content, "file content")
        self.mock_session.get.assert_called_once_with(
            "http://example.com/file.txt", timeout=self.client.timeout, stream=True
        )
        mock_response.raise_for_status.assert_called_once()

    def test_fetch_file_content_failure(self):
        mock_response = Mock()
        mock_response.status_code = 500
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError()
        self.mock_session.get.return_value = mock_response
        with self.assertRaises(requests.exceptions.HTTPError):
            self.client._fetch_file_content("http://example.com/file_error.txt")

    @patch.object(AFDBClient, "_make_api_request")
    def test_get_alphafold_entry_details_success(self, mock_api_request):
        mock_api_request.return_value = SAMPLE_API_RESPONSE_ENTRY
        details = self.client.get_alphafold_entry_details(SAMPLE_UNIPROT_ID)
        self.assertIsNotNone(details)
        self.assertEqual(details["uniprotAccession"], SAMPLE_UNIPROT_ID)
        self.assertEqual(details["cifUrl"], SAMPLE_API_RESPONSE_ENTRY[0]["cifUrl"])
        mock_api_request.assert_called_once_with(f"/prediction/{SAMPLE_UNIPROT_ID}")

    @patch.object(AFDBClient, "_make_api_request")
    def test_get_alphafold_entry_details_no_entry(self, mock_api_request):
        mock_api_request.return_value = []
        details = self.client.get_alphafold_entry_details("NONEXISTENT")
        self.assertIsNone(details)

    @patch.object(AFDBClient, "_make_api_request")
    def test_get_alphafold_entry_details_api_error(self, mock_api_request):
        mock_api_request.side_effect = requests.exceptions.RequestException("API down")
        details = self.client.get_alphafold_entry_details(SAMPLE_UNIPROT_ID)
        self.assertIsNone(details)

    def test_parse_plddt_from_cif_content_success(self):
        # Assumes method is named _parse_plddt_from_cif_content
        plddt_scores = self.client._parse_plddt_from_cif_content(
            SAMPLE_CIF_CONTENT, SAMPLE_UNIPROT_ID
        )
        expected_scores = [
            {"residue_number": 1, "plddt": 90.5},
            {"residue_number": 2, "plddt": 85.3},
            {"residue_number": 3, "plddt": 95.0},
        ]
        self.assertEqual(plddt_scores, expected_scores)

    def test_parse_plddt_from_cif_content_no_tags(self):
        cif_content_no_tags = "data_test\n#\nloop_\n_some_other_tag\n1\n"
        plddt_scores = self.client._parse_plddt_from_cif_content(
            cif_content_no_tags, SAMPLE_UNIPROT_ID
        )
        self.assertEqual(plddt_scores, [])

    def test_parse_plddt_from_cif_content_mismatched_lengths(self):
        plddt_scores = self.client._parse_plddt_from_cif_content(
            SAMPLE_CIF_CONTENT_MISMATCHED_LENGTHS, SAMPLE_UNIPROT_ID
        )
        self.assertEqual(plddt_scores, [])

    def test_parse_plddt_from_cif_content_bad_values(self):
        plddt_scores = self.client._parse_plddt_from_cif_content(
            SAMPLE_CIF_CONTENT_BAD_VALUES, SAMPLE_UNIPROT_ID
        )
        expected_scores = [{"residue_number": 3, "plddt": 95.0}]  # Only the valid one
        self.assertEqual(plddt_scores, expected_scores)

    def test_parse_plddt_from_cif_content_empty_input(self):
        plddt_scores = self.client._parse_plddt_from_cif_content("", SAMPLE_UNIPROT_ID)
        self.assertEqual(plddt_scores, [])

    @patch("gemmi.cif.read_string")
    def test_parse_plddt_from_cif_content_gemmi_error(self, mock_read_string):
        mock_read_string.side_effect = Exception("Gemmi parsing failed")
        plddt_scores = self.client._parse_plddt_from_cif_content(
            "corrupted cif", SAMPLE_UNIPROT_ID
        )
        self.assertEqual(plddt_scores, [])

    def test_parse_alphamissense_csv_content_success(self):
        # Assumes method is named _parse_alphamissense_csv_content
        am_scores = self.client._parse_alphamissense_csv_content(
            SAMPLE_ALPHAMISSENSE_CSV_CONTENT, SAMPLE_UNIPROT_ID
        )
        expected_scores = [
            {
                "residue_number": 1,
                "ref_aa": "M",
                "alt_aa": "A",
                "am_pathogenicity": 0.123,
                "am_class": "likely_benign",
            },
            {
                "residue_number": 2,
                "ref_aa": "L",
                "alt_aa": "P",
                "am_pathogenicity": 0.890,
                "am_class": "likely_pathogenic",
            },
            {
                "residue_number": 3,
                "ref_aa": "P",
                "alt_aa": "X",
                "am_pathogenicity": 0.500,
                "am_class": "ambiguous",
            },
        ]
        self.assertEqual(am_scores, expected_scores)

    def test_parse_alphamissense_csv_content_bad_header(self):
        am_scores = self.client._parse_alphamissense_csv_content(
            SAMPLE_ALPHAMISSENSE_CSV_CONTENT_BAD_HEADER, SAMPLE_UNIPROT_ID
        )
        self.assertEqual(am_scores, [])

    def test_parse_alphamissense_csv_content_short_row(self):
        am_scores = self.client._parse_alphamissense_csv_content(
            SAMPLE_ALPHAMISSENSE_CSV_CONTENT_SHORT_ROW, SAMPLE_UNIPROT_ID
        )
        self.assertEqual(am_scores, [])

    def test_parse_alphamissense_csv_content_invalid_variant_format(self):
        am_scores = self.client._parse_alphamissense_csv_content(
            SAMPLE_ALPHAMISSENSE_CSV_CONTENT_INVALID_VARIANT, SAMPLE_UNIPROT_ID
        )
        self.assertEqual(am_scores, [])

    def test_parse_alphamissense_csv_content_empty_after_header(self):
        csv_content = "protein_variant,am_pathogenicity,am_class\n"
        am_scores = self.client._parse_alphamissense_csv_content(
            csv_content, SAMPLE_UNIPROT_ID
        )
        self.assertEqual(am_scores, [])

    def test_parse_alphamissense_csv_content_empty_file(self):
        am_scores = self.client._parse_alphamissense_csv_content("", SAMPLE_UNIPROT_ID)
        self.assertEqual(am_scores, [])

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    @patch.object(AFDBClient, "_parse_plddt_from_cif_content")
    @patch.object(AFDBClient, "_parse_alphamissense_csv_content")
    def test_get_alphafold_data_plddt_success(
        self, mock_parse_am, mock_parse_plddt, mock_fetch, mock_get_details
    ):
        mock_get_details.return_value = {
            "cifUrl": "http://example.com/model.cif",
            "amAnnotationsUrl": None,
        }
        mock_fetch.return_value = "cif_file_content_mock"
        mock_parse_plddt.return_value = [{"residue_number": 1, "plddt": 90.0}]
        data = self.client.get_alphafold_data(SAMPLE_UNIPROT_ID, ["plddt"])
        self.assertIn("plddt", data)
        self.assertEqual(data["plddt"], [{"residue_number": 1, "plddt": 90.0}])
        mock_get_details.assert_called_once_with(SAMPLE_UNIPROT_ID)
        mock_fetch.assert_called_once_with("http://example.com/model.cif")
        mock_parse_plddt.assert_called_once_with(
            "cif_file_content_mock", SAMPLE_UNIPROT_ID
        )
        mock_parse_am.assert_not_called()

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    @patch.object(AFDBClient, "_parse_plddt_from_cif_content")
    @patch.object(AFDBClient, "_parse_alphamissense_csv_content")
    def test_get_alphafold_data_alphamissense_success(
        self, mock_parse_am, mock_parse_plddt, mock_fetch, mock_get_details
    ):
        mock_get_details.return_value = {
            "cifUrl": None,
            "amAnnotationsUrl": "http://example.com/am.csv.gz",
            "taxId": 9606,
        }
        mock_fetch.return_value = "am_file_content_mock"
        mock_parse_am.return_value = [{"residue_number": 1, "am_pathogenicity": 0.5}]
        data = self.client.get_alphafold_data(SAMPLE_UNIPROT_ID, ["alphamissense"])
        self.assertIn("alphamissense", data)
        self.assertEqual(
            data["alphamissense"], [{"residue_number": 1, "am_pathogenicity": 0.5}]
        )
        mock_get_details.assert_called_once_with(SAMPLE_UNIPROT_ID)
        mock_fetch.assert_called_once_with("http://example.com/am.csv.gz")
        mock_parse_am.assert_called_once_with("am_file_content_mock", SAMPLE_UNIPROT_ID)
        mock_parse_plddt.assert_not_called()

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    def test_get_alphafold_data_no_entry_details(self, mock_fetch, mock_get_details):
        mock_get_details.return_value = None
        data = self.client.get_alphafold_data(
            SAMPLE_UNIPROT_ID, ["plddt", "alphamissense"]
        )
        self.assertEqual(data, {"plddt": [], "alphamissense": []})
        mock_get_details.assert_called_once_with(SAMPLE_UNIPROT_ID)
        mock_fetch.assert_not_called()

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    def test_get_alphafold_data_plddt_no_cif_url(self, mock_fetch, mock_get_details):
        mock_get_details.return_value = {"cifUrl": None, "amAnnotationsUrl": None}
        data = self.client.get_alphafold_data(SAMPLE_UNIPROT_ID, ["plddt"])
        self.assertEqual(data["plddt"], [])
        mock_fetch.assert_not_called()

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    def test_get_alphafold_data_alphamissense_no_am_url_human(
        self, mock_fetch, mock_get_details
    ):
        mock_get_details.return_value = {
            "cifUrl": None,
            "amAnnotationsUrl": None,
            "taxId": 9606,
        }
        data = self.client.get_alphafold_data(SAMPLE_UNIPROT_ID, ["alphamissense"])
        self.assertEqual(data["alphamissense"], [])
        mock_fetch.assert_not_called()

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    def test_get_alphafold_data_alphamissense_no_am_url_non_human(
        self, mock_fetch, mock_get_details
    ):
        mock_get_details.return_value = {
            "cifUrl": None,
            "amAnnotationsUrl": None,
            "taxId": 10090,
        }
        data = self.client.get_alphafold_data(SAMPLE_UNIPROT_ID, ["alphamissense"])
        self.assertEqual(data["alphamissense"], [])
        mock_fetch.assert_not_called()

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    def test_get_alphafold_data_plddt_fetch_cif_fails(
        self, mock_fetch, mock_get_details
    ):
        mock_get_details.return_value = {"cifUrl": "http://example.com/model.cif"}
        mock_fetch.side_effect = requests.exceptions.RequestException("CIF Fetch error")
        data = self.client.get_alphafold_data(SAMPLE_UNIPROT_ID, ["plddt"])
        self.assertEqual(data["plddt"], [])
        mock_fetch.assert_called_once_with("http://example.com/model.cif")

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    def test_get_alphafold_data_am_fetch_fails(self, mock_fetch, mock_get_details):
        mock_get_details.return_value = {
            "amAnnotationsUrl": "http://example.com/am.csv.gz",
            "taxId": 9606,
        }
        mock_fetch.side_effect = requests.exceptions.RequestException("AM Fetch error")
        data = self.client.get_alphafold_data(SAMPLE_UNIPROT_ID, ["alphamissense"])
        self.assertEqual(data["alphamissense"], [])
        mock_fetch.assert_called_once_with("http://example.com/am.csv.gz")

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    @patch.object(AFDBClient, "_parse_plddt_from_cif_content")
    def test_get_alphafold_data_plddt_parse_fails(
        self, mock_parse_plddt, mock_fetch, mock_get_details
    ):
        mock_get_details.return_value = {"cifUrl": "http://example.com/model.cif"}
        mock_fetch.return_value = "cif_content_mock"
        mock_parse_plddt.return_value = []
        data = self.client.get_alphafold_data(SAMPLE_UNIPROT_ID, ["plddt"])
        self.assertEqual(data["plddt"], [])
        mock_parse_plddt.assert_called_once_with("cif_content_mock", SAMPLE_UNIPROT_ID)

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    @patch.object(AFDBClient, "_parse_alphamissense_csv_content")
    def test_get_alphafold_data_am_parse_fails(
        self, mock_parse_am, mock_fetch, mock_get_details
    ):
        mock_get_details.return_value = {
            "amAnnotationsUrl": "http://example.com/am.csv.gz",
            "taxId": 9606,
        }
        mock_fetch.return_value = "am_content_mock"
        mock_parse_am.return_value = []
        data = self.client.get_alphafold_data(SAMPLE_UNIPROT_ID, ["alphamissense"])
        self.assertEqual(data["alphamissense"], [])
        mock_parse_am.assert_called_once_with("am_content_mock", SAMPLE_UNIPROT_ID)

    @patch.object(AFDBClient, "get_alphafold_entry_details")
    @patch.object(AFDBClient, "_fetch_file_content")
    @patch.object(AFDBClient, "_parse_plddt_from_cif_content")
    @patch.object(AFDBClient, "_parse_alphamissense_csv_content")
    def test_get_alphafold_data_both_requested_one_fails_to_fetch(
        self, mock_parse_am, mock_parse_plddt, mock_fetch, mock_get_details
    ):
        mock_get_details.return_value = {
            "cifUrl": "http://example.com/model.cif",
            "amAnnotationsUrl": "http://example.com/am.csv.gz",
            "taxId": 9606,
        }

        def fetch_side_effect(url, **kwargs):
            if "model.cif" in url:
                raise requests.exceptions.RequestException("CIF Fetch error")
            elif "am.csv.gz" in url:
                return "am_content_mock"
            return None

        mock_fetch.side_effect = fetch_side_effect
        mock_parse_am.return_value = [{"residue_number": 1, "am_pathogenicity": 0.5}]

        data = self.client.get_alphafold_data(
            SAMPLE_UNIPROT_ID, ["plddt", "alphamissense"]
        )

        self.assertEqual(data["plddt"], [])
        self.assertEqual(
            data["alphamissense"], [{"residue_number": 1, "am_pathogenicity": 0.5}]
        )
        mock_get_details.assert_called_once_with(SAMPLE_UNIPROT_ID)
        self.assertIn(call("http://example.com/model.cif"), mock_fetch.call_args_list)
        self.assertIn(call("http://example.com/am.csv.gz"), mock_fetch.call_args_list)
        mock_parse_plddt.assert_not_called()
        mock_parse_am.assert_called_once_with("am_content_mock", SAMPLE_UNIPROT_ID)


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False)
