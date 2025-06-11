import logging
import unittest
from unittest.mock import Mock, patch

import requests

# Assuming the TEDClient is in this path for testing.
from protviz.data_retrieval.ted_client import TEDClient

# --- Sample data for mocking ---
SAMPLE_UNIPROT_ID = "Q53XC5"

# Sample successful response from the TED API
SAMPLE_TED_API_RESPONSE = {
    "data": [
        {
            "consensus_level": "High",
            "chopping": "1-124_CATH-3.30.240.10",
            "num_segments": 1,
            "cath_label": "CATH-3.30.240.10",
            "nres_domain": 124,
        },
        {
            "consensus_level": "Medium",
            "chopping": "130-250_CATH-1.10.510.10",
            "num_segments": 1,
            "cath_label": "CATH-1.10.510.10",
            "nres_domain": 121,
        },
    ]
}

# Sample response with no 'data' key
EMPTY_TED_API_RESPONSE = {"message": "No data found"}

# Sample response with 'data' key but not a list
MALFORMED_TED_API_RESPONSE = {"data": {"error": "Unexpected format"}}


class TestTEDClient(unittest.TestCase):
    """Unit tests for the TEDClient class."""

    def setUp(self):
        """Set up a test client with patched dependencies."""
        client_module_path = "protviz.data_retrieval.ted_client"

        # Patch dependencies to prevent side effects
        self.platformdirs_patcher = patch(
            f"{client_module_path}.platformdirs.user_cache_dir",
            return_value="/fake/cache/dir",
        )
        self.makedirs_patcher = patch(f"{client_module_path}.os.makedirs")
        self.requests_cache_patcher = patch(
            f"{client_module_path}.requests_cache.CachedSession"
        )

        # Start patchers
        self.mock_user_cache_dir = self.platformdirs_patcher.start()
        self.mock_makedirs = self.makedirs_patcher.start()
        self.mock_cached_session_cls = self.requests_cache_patcher.start()

        # Instantiate the client; it will use the mocked session
        self.client = TEDClient(timeout=5)
        self.mock_session = self.client.session

        # Suppress logging for clean test output
        logging.disable(logging.CRITICAL)

    def tearDown(self):
        """Stop all patchers and re-enable logging."""
        self.platformdirs_patcher.stop()
        self.makedirs_patcher.stop()
        self.requests_cache_patcher.stop()
        logging.disable(logging.NOTSET)

    def test_make_request_success(self):
        """Test a successful API request."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = SAMPLE_TED_API_RESPONSE
        mock_response.content = b'{"data": []}'  # Ensure content is not empty
        self.mock_session.get.return_value = mock_response

        response = self.client._make_request("summary", SAMPLE_UNIPROT_ID)

        self.assertEqual(response, SAMPLE_TED_API_RESPONSE)
        expected_url = f"{self.client.TED_API_URL}summary/{SAMPLE_UNIPROT_ID}"
        self.mock_session.get.assert_called_once_with(expected_url, timeout=5)
        mock_response.raise_for_status.assert_called_once()

    def test_make_request_http_error(self):
        """Test that an HTTPError is raised for bad status codes."""
        mock_response = Mock()
        mock_response.status_code = 404
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
            response=mock_response
        )
        self.mock_session.get.return_value = mock_response

        with self.assertRaises(requests.exceptions.HTTPError):
            self.client._make_request("summary", "NOT_FOUND_ID")

    def test_make_request_json_decode_error(self):
        """Test handling of invalid JSON responses."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b"not a valid json"
        mock_response.text = "not a valid json"
        mock_response.json.side_effect = requests.exceptions.JSONDecodeError(
            "msg", "doc", 0
        )
        self.mock_session.get.return_value = mock_response

        with self.assertRaises(ValueError):
            self.client._make_request("summary", "BAD_JSON_ID")

    def test_make_request_request_exception(self):
        """Test handling of general request exceptions."""
        self.mock_session.get.side_effect = requests.exceptions.RequestException(
            "Network failure"
        )

        with self.assertRaises(requests.exceptions.RequestException):
            self.client._make_request("summary", "NETWORK_ERROR_ID")

    def test_make_request_empty_content(self):
        """Test that an empty response content is handled correctly."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b""  # Empty content
        self.mock_session.get.return_value = mock_response

        response = self.client._make_request("summary", "EMPTY_CONTENT_ID")
        self.assertEqual(response, {})

    @patch.object(TEDClient, "_make_request")
    def test_get_ted_annotations_success(self, mock_make_request):
        """Test successful parsing of TED annotations."""
        mock_make_request.return_value = SAMPLE_TED_API_RESPONSE

        annotations = self.client.get_TED_annotations(SAMPLE_UNIPROT_ID)

        self.assertEqual(len(annotations), 2)

        # Check the first annotation
        self.assertEqual(annotations[0]["uniprot_acc"], SAMPLE_UNIPROT_ID)
        self.assertEqual(annotations[0]["consensus_level"], "High")
        self.assertEqual(annotations[0]["chopping"], "1-124_CATH-3.30.240.10")
        self.assertEqual(annotations[0]["cath_label"], "CATH-3.30.240.10")

        # Check the second annotation
        self.assertEqual(annotations[1]["num_segments"], 1)
        self.assertEqual(annotations[1]["nres_domain"], 121)

        mock_make_request.assert_called_once_with("summary", SAMPLE_UNIPROT_ID)

    @patch.object(TEDClient, "_make_request")
    def test_get_ted_annotations_empty_response(self, mock_make_request):
        """Test handling of an API response with no 'data' key."""
        mock_make_request.return_value = EMPTY_TED_API_RESPONSE
        annotations = self.client.get_TED_annotations(SAMPLE_UNIPROT_ID)
        self.assertEqual(annotations, [])

    @patch.object(TEDClient, "_make_request")
    def test_get_ted_annotations_malformed_response(self, mock_make_request):
        """Test handling of a response where 'data' is not a list."""
        mock_make_request.return_value = MALFORMED_TED_API_RESPONSE
        annotations = self.client.get_TED_annotations(SAMPLE_UNIPROT_ID)
        self.assertEqual(annotations, [])

    @patch.object(TEDClient, "_make_request")
    def test_get_ted_annotations_request_fails(self, mock_make_request):
        """Test that an empty list is returned if the request fails."""
        mock_make_request.side_effect = requests.exceptions.RequestException(
            "API is down"
        )
        annotations = self.client.get_TED_annotations(SAMPLE_UNIPROT_ID)
        self.assertEqual(annotations, [])

    @patch.object(TEDClient, "_make_request")
    def test_get_ted_annotations_no_valid_items_in_data(self, mock_make_request):
        """Test a response where 'data' contains no valid dictionary items."""
        mock_make_request.return_value = {"data": ["just a string", None, 123]}
        annotations = self.client.get_TED_annotations(SAMPLE_UNIPROT_ID)
        self.assertEqual(annotations, [])


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False)
