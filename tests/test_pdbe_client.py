import logging
import unittest
from unittest.mock import Mock, patch

import requests

# Assuming the PDBeClient is in this path for testing.
from protviz.data_retrieval.pdbe_client import PDBeClient

# --- Sample data for mocking ---
SAMPLE_UNIPROT_ID = "P00533"  # EGFR, a good example with ligands and structures

# Sample response for get_pdb_coverage
SAMPLE_COVERAGE_RESPONSE = {
    SAMPLE_UNIPROT_ID: [
        {"pdb_id": "2gs2", "unp_start": 696, "unp_end": 1022},
        {"pdb_id": "2gs7", "unp_start": 696, "unp_end": 1022},
        # Entry with missing keys to test robustness
        {"unp_start": 1, "unp_end": 100},
    ]
}

# Sample response for get_pdb_ligand_interactions
SAMPLE_LIGAND_RESPONSE = {
    SAMPLE_UNIPROT_ID: {
        "dataType": "LIGAND BINDING SITES",
        "data": [
            {
                "accession": "ATP",
                "pdb_id": "2gs2",
                "residues": [
                    {
                        "startIndex": 721,
                        "endIndex": 721,
                        "startCode": "L",
                        "endCode": "L",
                    },
                    {
                        "startIndex": 793,
                        "endIndex": 793,
                        "startCode": "M",
                        "endCode": "M",
                    },
                ],
            },
            {
                "accession": "MG",
                "pdb_id": "2gs2",
                "residues": [
                    {
                        "startIndex": 795,
                        "endIndex": 795,
                        "startCode": "D",
                        "endCode": "D",
                    }
                ],
            },
            # Malformed residue entry to test robustness
            {"accession": "HOH", "pdb_id": "2gs2", "residues": [{"startCode": "X"}]},
        ],
    }
}


class TestPDBeClient(unittest.TestCase):
    """Unit tests for the PDBeClient class."""

    def setUp(self):
        """Set up a test client with patched dependencies."""
        client_module_path = "protviz.data_retrieval.pdbe_client"

        self.platformdirs_patcher = patch(
            f"{client_module_path}.platformdirs.user_cache_dir",
            return_value="/fake/cache/dir",
        )
        self.makedirs_patcher = patch(f"{client_module_path}.os.makedirs")
        self.requests_cache_patcher = patch(
            f"{client_module_path}.requests_cache.CachedSession"
        )

        self.mock_user_cache_dir = self.platformdirs_patcher.start()
        self.mock_makedirs = self.makedirs_patcher.start()
        self.mock_cached_session_cls = self.requests_cache_patcher.start()

        self.client = PDBeClient(timeout=5)
        self.mock_session = self.client.session

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
        mock_response.json.return_value = {"data": "success"}
        mock_response.content = b'{"data": "success"}'
        self.mock_session.get.return_value = mock_response

        response = self.client._make_request("/test", SAMPLE_UNIPROT_ID)

        self.assertEqual(response, {"data": "success"})
        expected_url = f"{self.client.PDBE_GRAPH_API_BASE_URL}/test/{SAMPLE_UNIPROT_ID}"
        self.mock_session.get.assert_called_once_with(expected_url, timeout=5)
        mock_response.raise_for_status.assert_called_once()

    def test_make_request_http_error(self):
        """Test that HTTPError is raised for bad status codes."""
        mock_response = Mock()
        mock_response.status_code = 500
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
            response=mock_response
        )
        self.mock_session.get.return_value = mock_response

        with self.assertRaises(requests.exceptions.HTTPError):
            self.client._make_request("/error", "ERROR_ID")

    def test_make_request_json_decode_error(self):
        """Test handling of invalid JSON responses."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b"not valid json"
        mock_response.text = "not valid json"
        mock_response.json.side_effect = requests.exceptions.JSONDecodeError(
            "msg", "doc", 0
        )
        self.mock_session.get.return_value = mock_response

        with self.assertRaises(ValueError):
            self.client._make_request("/badjson", "BAD_JSON_ID")

    @patch.object(PDBeClient, "_make_request")
    def test_get_pdb_coverage_success(self, mock_make_request):
        """Test successful fetching and parsing of PDB coverage."""
        mock_make_request.return_value = SAMPLE_COVERAGE_RESPONSE

        coverage = self.client.get_pdb_coverage(SAMPLE_UNIPROT_ID)

        self.assertEqual(len(coverage), 2)  # Should skip the malformed entry
        self.assertEqual(coverage[0]["pdb_id"], "2gs2")
        self.assertEqual(coverage[0]["unp_start"], 696)
        self.assertEqual(coverage[1]["unp_end"], 1022)
        mock_make_request.assert_called_once_with(
            "/mappings/best_structures", SAMPLE_UNIPROT_ID
        )

    @patch.object(PDBeClient, "_make_request")
    def test_get_pdb_coverage_request_fails(self, mock_make_request):
        """Test that an empty list is returned if the request fails."""
        mock_make_request.side_effect = requests.exceptions.RequestException
        coverage = self.client.get_pdb_coverage(SAMPLE_UNIPROT_ID)
        self.assertEqual(coverage, [])

    @patch.object(PDBeClient, "_make_request")
    def test_get_pdb_coverage_no_data(self, mock_make_request):
        """Test response with no data for the given UniProt ID."""
        mock_make_request.return_value = {"OTHER_ID": []}
        coverage = self.client.get_pdb_coverage(SAMPLE_UNIPROT_ID)
        self.assertEqual(coverage, [])

    @patch.object(PDBeClient, "_make_request")
    def test_get_pdb_ligand_interactions_success(self, mock_make_request):
        """Test successful fetching and parsing of ligand interactions."""
        mock_make_request.return_value = SAMPLE_LIGAND_RESPONSE

        interactions = self.client.get_pdb_ligand_interactions(SAMPLE_UNIPROT_ID)

        self.assertEqual(len(interactions), 2)  # ATP and MG ligands

        atp_interaction = next(i for i in interactions if i["ligand_id"] == "ATP")
        self.assertEqual(len(atp_interaction["binding_site_uniprot_residues"]), 2)
        self.assertEqual(
            atp_interaction["binding_site_uniprot_residues"][0]["startIndex"], 721
        )
        self.assertEqual(atp_interaction["pdb_id"], "2gs2")

        mg_interaction = next(i for i in interactions if i["ligand_id"] == "MG")
        self.assertEqual(len(mg_interaction["binding_site_uniprot_residues"]), 1)

        mock_make_request.assert_called_once_with(
            "/uniprot/ligand_sites", SAMPLE_UNIPROT_ID
        )

    @patch.object(PDBeClient, "_make_request")
    def test_get_pdb_ligand_interactions_request_fails(self, mock_make_request):
        """Test that an empty list is returned if ligand request fails."""
        mock_make_request.side_effect = requests.exceptions.RequestException
        interactions = self.client.get_pdb_ligand_interactions(SAMPLE_UNIPROT_ID)
        self.assertEqual(interactions, [])

    @patch.object(PDBeClient, "_make_request")
    def test_get_pdb_ligand_interactions_malformed_response(self, mock_make_request):
        """Test various forms of malformed ligand data."""
        # Case 1: Wrong data type
        mock_make_request.return_value = {
            SAMPLE_UNIPROT_ID: {"dataType": "WRONG_TYPE", "data": []}
        }
        self.assertEqual(self.client.get_pdb_ligand_interactions(SAMPLE_UNIPROT_ID), [])

        # Case 2: 'data' is not a list
        mock_make_request.return_value = {
            SAMPLE_UNIPROT_ID: {
                "dataType": "LIGAND BINDING SITES",
                "data": {"not": "a list"},
            }
        }
        self.assertEqual(self.client.get_pdb_ligand_interactions(SAMPLE_UNIPROT_ID), [])

        # Case 3: 'residues' is not a list
        mock_make_request.return_value = {
            SAMPLE_UNIPROT_ID: {
                "dataType": "LIGAND BINDING SITES",
                "data": [{"accession": "LIG", "residues": "not a list"}],
            }
        }
        self.assertEqual(self.client.get_pdb_ligand_interactions(SAMPLE_UNIPROT_ID), [])

    def test_clear_cache(self):
        """Test that the cache clear method is called on the session."""
        self.client.clear_cache()
        self.mock_session.cache.clear.assert_called_once()


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False)
