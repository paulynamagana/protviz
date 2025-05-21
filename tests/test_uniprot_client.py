from unittest.mock import MagicMock

import pytest
import requests

from protviz.data_retrieval.uniprot_client import get_protein_sequence_length


# -- Test cases for get_protein_sequence_length function --
def test_get_protein_sequence_length_success(mocker):
    """
    Test get_protein_sequence_length with a valid UniProt ID.
    It should return the correct sequence length.
    """

    uniprot_id = "P0DTC2"
    expected_length = 1273
    mock_api_response_json = {
        "sequence": {
            "length": expected_length,
            "value": "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHV",
        }
        # add any other necessary fields here!
    }
    # Create a mock response object for requests.get()
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = mock_api_response_json
    mock_response.raise_for_status.return_value = None  # Simulate no HTTP error

    # Patch 'requests.get' within the scope of your uniprot_client module
    # This tells pytest to replace requests.get with our mock_response when it's called
    # from within protviz.data_retrieval.uniprot_client
    mocker.patch(
        "protviz.data_retrieval.uniprot_client.requests.get", return_value=mock_response
    )

    # Call the function under test
    actual_length = get_protein_sequence_length(uniprot_id)

    # Assertions: Check if the function behaved as expected
    assert actual_length == expected_length

    # Example of checking call:
    # protviz.data_retrieval.uniprot_client.requests.get.assert_called_once_with(
    #     f"https://rest.uniprot.org/uniprotkb/{uniprot_id}",
    #     headers={"Accept": "application/json"}

    def test_get_protein_sequence_length_not_found(mocker):
        """
        Test get_protein_sequence_length when the UniProt ID is not found (404 error).
        It should raise a ValueError.
        """
        uniprot_id = "NOTANID"

        # Create a mock response object for a 404 error
        mock_response = MagicMock()
        mock_response.status_code = 404
        # Simulate raise_for_status() raising an HTTPError for 404
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
            "404 Client Error: Not Found for url", response=mock_response
        )

        mocker.patch(
            "protviz.data_retrieval.uniprot_client.requests.get",
            return_value=mock_response,
        )

        # Use pytest.raises to assert that a specific exception is raised
        with pytest.raises(ValueError) as excinfo:
            get_protein_sequence_length(uniprot_id)

        # Optionally, check the exception message
        assert f"Protein with UniProt ID '{uniprot_id}' not found." in str(
            excinfo.value
        )

    def test_get_protein_sequence_length_api_error(mocker):
        """
        Test get_protein_sequence_length for a generic API error (e.g., 500).
        It should re-raise the requests.exceptions.HTTPError.
        """
        uniprot_id = "P12345"

        mock_response = MagicMock()
        mock_response.status_code = 500
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
            "500 Server Error", response=mock_response
        )

        mocker.patch(
            "protviz.data_retrieval.uniprot_client.requests.get",
            return_value=mock_response,
        )

        with pytest.raises(requests.exceptions.HTTPError) as excinfo:
            get_protein_sequence_length(uniprot_id)

        assert "500 Server Error" in str(excinfo.value)

    def test_get_protein_sequence_length_missing_sequence_key(mocker):
        """
        Test get_protein_sequence_length when the API response is valid JSON
        but missing the expected 'sequence' or 'length' key.
        It should raise a KeyError.
        """
        uniprot_id = "P67890"
        mock_api_response_json = {
            "accession": uniprot_id,
            "id": "SOME_ID_HUMAN",
            # Missing 'sequence' key
        }

        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = mock_api_response_json
        mock_response.raise_for_status.return_value = None

        mocker.patch(
            "protviz.data_retrieval.uniprot_client.requests.get",
            return_value=mock_response,
        )

        with pytest.raises(KeyError) as excinfo:
            get_protein_sequence_length(uniprot_id)

        assert "Sequence information not found in the response." in str(excinfo.value)


# To run these tests:
# Run: pytest
#
# Check a scpecific test file:
# pytest tests/data_retrieval/test_uniprot_client.py
#
# Or a specific test function:
# pytest tests/data_retrieval/test_uniprot_client.py::test_get_protein_sequence_length_success
