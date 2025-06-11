import unittest

import requests
import requests_mock

from protviz.data_retrieval.uniprot_client import get_protein_sequence_length


class TestUniprotClient(unittest.TestCase):
    """Unit tests for the UniProt client."""

    UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb"

    @requests_mock.Mocker()
    def test_get_protein_sequence_length_success(self, m):
        """Test a successful retrieval of protein sequence length."""
        uniprot_id = "P02768"
        mock_url = f"{self.UNIPROT_API_URL}/{uniprot_id}"
        mock_response = {"sequence": {"length": 609}}

        m.get(mock_url, json=mock_response, status_code=200)

        length = get_protein_sequence_length(uniprot_id)
        self.assertEqual(length, 609)
        self.assertTrue(m.called)

    @requests_mock.Mocker()
    def test_get_protein_sequence_length_not_found(self, m):
        """Test the function's handling of a 404 Not Found error."""
        uniprot_id = "INVALIDID"
        mock_url = f"{self.UNIPROT_API_URL}/{uniprot_id}"

        m.get(mock_url, status_code=404)

        with self.assertRaises(ValueError) as context:
            get_protein_sequence_length(uniprot_id)

        self.assertIn(f"UniProt ID '{uniprot_id}' not found.", str(context.exception))

    @requests_mock.Mocker()
    def test_get_protein_sequence_length_server_error(self, m):
        """Test the function's handling of a 500 server error."""
        uniprot_id = "P0DTC2"
        mock_url = f"{self.UNIPROT_API_URL}/{uniprot_id}"

        m.get(mock_url, status_code=500)

        with self.assertRaises(requests.exceptions.HTTPError):
            get_protein_sequence_length(uniprot_id)

    @requests_mock.Mocker()
    def test_get_protein_sequence_length_malformed_response(self, m):
        """Test the function's handling of a response missing sequence data."""
        uniprot_id = "P12345"
        mock_url = f"{self.UNIPROT_API_URL}/{uniprot_id}"
        mock_response = {"entry": "some other data"}  # No 'sequence' key

        m.get(mock_url, json=mock_response, status_code=200)

        with self.assertRaises(KeyError) as context:
            get_protein_sequence_length(uniprot_id)

        self.assertIn("Sequence information not found", str(context.exception))


if __name__ == "__main__":
    unittest.main()
