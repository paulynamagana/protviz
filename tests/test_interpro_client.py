import logging
import unittest
from unittest.mock import Mock, patch

import requests

from protviz.data_retrieval.interpro_client import InterProClient

# Sample data for mocking
SAMPLE_UNIPROT_ID = "P04637"
SAMPLE_UNIPROT_ID_LOWER = (
    "p04637"  # Test case insensitivity for protein accession matching
)

# A more complete and structured sample response from InterPro API
# /api/entry/interpro/protein/uniprot/{accession}
SAMPLE_INTERPRO_SUMMARY_RESPONSE = {
    "count": 1,
    "next": None,
    "previous": None,
    "results": [
        {
            "metadata": {
                "accession": "IPR002117",
                "name": "P53 tumour suppressor family",
                "type": "Family",
                "source_database": "interpro",
                "member_databases": {
                    "pfam": {"PF00870": "P53 DNA-binding domain"},
                    "prints": {"PR00098": "P53"},
                    "cathgene3d": {"G3DSA:2.60.40.1230": "P53 DNA-binding domain-like"},
                },
                "go_terms": None,
                "description": None,
                "integrated": "IPR002117",  # Often the same as accession for InterPro entries
                "hierarchy": None,
            },
            "proteins": [
                {
                    "accession": SAMPLE_UNIPROT_ID_LOWER,  # API returns lowercase
                    "protein_length": 393,
                    "source_database": "uniprot",
                    "entry_protein_locations": [
                        {
                            "fragments": [
                                {"start": 102, "end": 292, "dc-status": "CONTINUOUS"},
                                {
                                    "start": 1,
                                    "end": 50,
                                    "dc-status": "CONTINUOUS",
                                },  # Another fragment for testing
                            ],
                            "model_acc": "PF00870",  # Relates to a member db signature, not InterPro entry itself
                        }
                    ],
                    "ida_accession": None,
                    "ida_id": None,
                    "entry_type": "Family",  # Redundant from metadata
                }
            ],
        },
        {  # Another InterPro entry to test multiple results and different member DBs
            "metadata": {
                "accession": "IPR011615",
                "name": "P53, tetramerisation motif",
                "type": "Domain",
                "source_database": "interpro",
                "member_databases": {
                    "pfam": {"PF07710": "P53 tetramerisation motif"},
                    "smart": {"SM00106": "P53T"},
                },
                "go_terms": None,
                "description": None,
                "integrated": "IPR011615",
                "hierarchy": None,
            },
            "proteins": [
                {
                    "accession": SAMPLE_UNIPROT_ID_LOWER,
                    "protein_length": 393,
                    "source_database": "uniprot",
                    "entry_protein_locations": [
                        {
                            "fragments": [
                                {"start": 325, "end": 356, "dc-status": "CONTINUOUS"}
                            ],
                            "model_acc": "PF07710",
                        }
                    ],
                }
            ],
        },
        {  # InterPro entry with NO matching member db for a specific test
            "metadata": {
                "accession": "IPR000000",
                "name": "No Target DB Entry",
                "type": "Family",
                "source_database": "interpro",
                "member_databases": {"otherdb": {"OTH123": "Other DB entry"}},
            },
            "proteins": [
                {
                    "accession": SAMPLE_UNIPROT_ID_LOWER,
                    "protein_length": 393,
                    "source_database": "uniprot",
                    "entry_protein_locations": [
                        {"fragments": [{"start": 10, "end": 20}]}
                    ],
                }
            ],
        },
        {  # InterPro entry with NO locations for the specific protein
            "metadata": {
                "accession": "IPR999999",
                "name": "No Locations Entry",
                "type": "Domain",
                "source_database": "interpro",
                "member_databases": {"pfam": {"PF99999": "No Location Pfam"}},
            },
            "proteins": [
                {  # Belongs to a different protein
                    "accession": "q12345",
                    "protein_length": 100,
                    "source_database": "uniprot",
                    "entry_protein_locations": [
                        {"fragments": [{"start": 1, "end": 10}]}
                    ],
                }
                # No entry for SAMPLE_UNIPROT_ID here
            ],
        },
    ],
}

MALFORMED_INTERPRO_RESPONSE = {
    "count": 1,
    "results": [
        {"metadata": "not_a_dict"}  # Malformed metadata
    ],
}


class TestInterProClient(unittest.TestCase):
    def setUp(self):
        client_module_path = "protviz.data_retrieval.interpro_client"

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

        self.client = InterProClient(timeout=5)
        self.mock_session = self.client.session

        logging.disable(logging.CRITICAL)

    def tearDown(self):
        self.platformdirs_patcher.stop()
        self.makedirs_patcher.stop()
        self.requests_cache_patcher.stop()
        logging.disable(logging.NOTSET)

    def test_fetch_protein_interpro_summary_success(self):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"data": "success"}
        mock_response.content = b'{"data": "success"}'
        self.mock_session.get.return_value = mock_response

        response = self.client._fetch_protein_interpro_summary(SAMPLE_UNIPROT_ID)
        self.assertEqual(response, {"data": "success"})

        expected_url = (
            f"{self.client.INTERPRO_API_BASE_URL}/protein/uniprot/{SAMPLE_UNIPROT_ID}"
        )
        self.mock_session.get.assert_called_once_with(
            expected_url,
            timeout=self.client.timeout,
            headers={"Accept": "application/json"},
        )
        mock_response.raise_for_status.assert_called_once()

    def test_fetch_protein_interpro_summary_http_error_404(self):
        mock_response = Mock()
        mock_response.status_code = 404
        http_error = requests.exceptions.HTTPError(response=mock_response)
        mock_response.raise_for_status.side_effect = http_error
        self.mock_session.get.return_value = mock_response

        response = self.client._fetch_protein_interpro_summary("NOT_FOUND_ID")
        self.assertIsNone(response)

    def test_fetch_protein_interpro_summary_http_error_other(self):
        mock_response = Mock()
        mock_response.status_code = 500
        http_error = requests.exceptions.HTTPError(response=mock_response)
        mock_response.raise_for_status.side_effect = http_error
        self.mock_session.get.return_value = mock_response

        with self.assertRaises(requests.exceptions.HTTPError):
            self.client._fetch_protein_interpro_summary("ERROR_ID")

    def test_fetch_protein_interpro_summary_json_decode_error(self):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b"not a json"
        mock_response.text = "not a json"
        mock_response.json.side_effect = requests.exceptions.JSONDecodeError(
            "msg", "doc", 0
        )
        self.mock_session.get.return_value = mock_response

        with self.assertRaises(ValueError):
            self.client._fetch_protein_interpro_summary("INVALID_JSON_ID")

    def test_fetch_protein_interpro_summary_request_exception(self):
        self.mock_session.get.side_effect = requests.exceptions.RequestException(
            "Network error"
        )
        with self.assertRaises(requests.exceptions.RequestException):
            self.client._fetch_protein_interpro_summary("NETWORK_ERROR_ID")

    def test_fetch_protein_interpro_summary_empty_response_content(self):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.content = b""
        self.mock_session.get.return_value = mock_response

        response = self.client._fetch_protein_interpro_summary("EMPTY_CONTENT_ID")
        self.assertIsNone(response)

    def test_extract_member_db_annotations_pfam_success(self):
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID, SAMPLE_INTERPRO_SUMMARY_RESPONSE, "pfam"
        )
        # CORRECTED: The sample response contains 3 InterPro entries that have a "pfam" member db.
        # Two have locations for P04637, and one does not, but is still included.
        self.assertEqual(len(annotations), 3)

        # Check entry 1 (IPR002117 -> PF00870)
        entry1 = next(ann for ann in annotations if ann["accession"] == "PF00870")
        self.assertEqual(entry1["interpro_accession"], "IPR002117")
        self.assertEqual(entry1["name"], "P53 tumour suppressor family")
        self.assertEqual(entry1["description"], "P53 DNA-binding domain")
        self.assertEqual(entry1["entry_type"], "Family")
        self.assertIn({"start": 102, "end": 292}, entry1["locations"])
        self.assertIn({"start": 1, "end": 50}, entry1["locations"])

        # Check entry 2 (IPR011615 -> PF07710)
        entry2 = next(ann for ann in annotations if ann["accession"] == "PF07710")
        self.assertEqual(entry2["interpro_accession"], "IPR011615")
        self.assertEqual(entry2["name"], "P53, tetramerisation motif")
        self.assertEqual(entry2["description"], "P53 tetramerisation motif")
        self.assertEqual(entry2["entry_type"], "Domain")
        self.assertEqual(entry2["locations"], [{"start": 325, "end": 356}])

        # Check entry 3 (IPR999999 -> PF99999)
        entry3 = next(ann for ann in annotations if ann["accession"] == "PF99999")
        self.assertEqual(entry3["interpro_accession"], "IPR999999")
        self.assertEqual(entry3["name"], "No Locations Entry")
        self.assertEqual(entry3["locations"], [])  # Crucially, locations are empty

    def test_extract_member_db_annotations_cathgene3d_success(self):
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID, SAMPLE_INTERPRO_SUMMARY_RESPONSE, "cathgene3d"
        )
        self.assertEqual(len(annotations), 1)
        entry1 = annotations[0]
        self.assertEqual(entry1["accession"], "G3DSA:2.60.40.1230")
        self.assertEqual(entry1["interpro_accession"], "IPR002117")
        self.assertIn({"start": 102, "end": 292}, entry1["locations"])

    def test_extract_member_db_annotations_no_summary_data(self):
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID, None, "pfam"
        )
        self.assertEqual(annotations, [])

    def test_extract_member_db_annotations_empty_results(self):
        empty_summary = {"results": []}
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID, empty_summary, "pfam"
        )
        self.assertEqual(annotations, [])

    def test_extract_member_db_annotations_no_results_key(self):
        no_results_summary = {"something_else": []}
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID, no_results_summary, "pfam"
        )
        self.assertEqual(annotations, [])

    def test_extract_member_db_annotations_malformed_entry_in_results(self):
        summary_with_malformed_entry = {"results": ["not_a_dict"]}
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID, summary_with_malformed_entry, "pfam"
        )
        self.assertEqual(annotations, [])

    def test_extract_member_db_annotations_malformed_metadata(self):
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID, MALFORMED_INTERPRO_RESPONSE, "pfam"
        )
        self.assertEqual(annotations, [])

    def test_extract_member_db_annotations_target_db_not_in_entry(self):
        single_entry_summary = {
            "results": [
                r
                for r in SAMPLE_INTERPRO_SUMMARY_RESPONSE["results"]
                if r["metadata"]["accession"] == "IPR000000"
            ]
        }
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID,
            single_entry_summary,
            "pfam",
        )
        self.assertEqual(annotations, [])

    def test_extract_member_db_annotations_no_locations_for_protein(self):
        # Using IPR999999 which has PFAM but no locations for SAMPLE_UNIPROT_ID
        summary_no_loc = {
            "results": [
                r
                for r in SAMPLE_INTERPRO_SUMMARY_RESPONSE["results"]
                if r["metadata"]["accession"] == "IPR999999"
            ]
        }
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID, summary_no_loc, "pfam"
        )
        # CORRECTED: The code creates an entry for the member signature (PF99999)
        # but correctly assigns an empty list for locations because SAMPLE_UNIPROT_ID
        # was not found in the 'proteins' list for this entry.
        self.assertEqual(len(annotations), 1)
        self.assertEqual(annotations[0]["accession"], "PF99999")
        self.assertEqual(annotations[0]["locations"], [])

    def test_extract_member_db_annotations_locations_with_bad_start_end(self):
        bad_loc_summary = {
            "results": [
                {
                    "metadata": {
                        "accession": "IPR002117",
                        "type": "Family",
                        "name": "Test",
                        "member_databases": {"pfam": {"PF00870": "Test Pfam"}},
                    },
                    "proteins": [
                        {
                            "accession": SAMPLE_UNIPROT_ID_LOWER,
                            "entry_protein_locations": [
                                {
                                    "fragments": [
                                        {"start": "abc", "end": 292},
                                        {"start": 102, "end": "xyz"},
                                        {"start": 200, "end": 300},
                                    ]
                                }
                            ],
                        }
                    ],
                }
            ]
        }
        annotations = self.client._extract_member_db_annotations(
            SAMPLE_UNIPROT_ID, bad_loc_summary, "pfam"
        )
        self.assertEqual(len(annotations), 1)
        self.assertEqual(annotations[0]["accession"], "PF00870")
        self.assertEqual(annotations[0]["locations"], [{"start": 200, "end": 300}])

    @patch.object(InterProClient, "_fetch_protein_interpro_summary")
    @patch.object(InterProClient, "_extract_member_db_annotations")
    def test_get_pfam_annotations_success(self, mock_extract, mock_fetch):
        mock_fetch.return_value = {"results": ["dummy_summary_data"]}
        mock_extract.return_value = [{"accession": "PF00870", "name": "Test Pfam"}]

        result = self.client.get_pfam_annotations(SAMPLE_UNIPROT_ID)

        mock_fetch.assert_called_once_with(SAMPLE_UNIPROT_ID)
        mock_extract.assert_called_once_with(
            SAMPLE_UNIPROT_ID, {"results": ["dummy_summary_data"]}, "pfam"
        )
        self.assertEqual(result, [{"accession": "PF00870", "name": "Test Pfam"}])

    @patch.object(InterProClient, "_fetch_protein_interpro_summary")
    def test_get_pfam_annotations_fetch_returns_none(self, mock_fetch):
        mock_fetch.return_value = None
        result = self.client.get_pfam_annotations(SAMPLE_UNIPROT_ID)
        self.assertEqual(result, [])
        mock_fetch.assert_called_once_with(SAMPLE_UNIPROT_ID)

    @patch.object(InterProClient, "_fetch_protein_interpro_summary")
    def test_get_pfam_annotations_fetch_raises_exception(self, mock_fetch):
        mock_fetch.side_effect = requests.exceptions.RequestException("API down")
        result = self.client.get_pfam_annotations(SAMPLE_UNIPROT_ID)
        self.assertEqual(result, [])
        mock_fetch.assert_called_once_with(SAMPLE_UNIPROT_ID)

    @patch.object(InterProClient, "_fetch_protein_interpro_summary")
    @patch.object(InterProClient, "_extract_member_db_annotations")
    def test_get_cathgene3d_annotations_success(self, mock_extract, mock_fetch):
        mock_fetch.return_value = {"results": ["dummy_summary_data_cath"]}
        mock_extract.return_value = [{"accession": "G3DSA:1.2.3", "name": "Test CATH"}]

        result = self.client.get_cathgene3d_annotations(SAMPLE_UNIPROT_ID)

        mock_fetch.assert_called_once_with(SAMPLE_UNIPROT_ID)
        mock_extract.assert_called_once_with(
            SAMPLE_UNIPROT_ID, {"results": ["dummy_summary_data_cath"]}, "cathgene3d"
        )
        self.assertEqual(result, [{"accession": "G3DSA:1.2.3", "name": "Test CATH"}])


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False)
