# src/protviz/data_retrieval/ted_client.py
import logging
from typing import Any, Dict, List

import requests

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class TEDClient:
    """
    A class to interact with the TED annotations
    """

    TED_API_URL = "https://ted.cathdb.info/api/v1/uniprot/"

    def __init__(self, timeout: int = 10):
        """
        Initialize the TED client with a timeout for requests.

        :param timeout: Timeout for requests in seconds
        """
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update(
            {"Accept": "application/json", "Content-Type": "application/json"}
        )

    def _make_request(self, endpoint: str, uniprot_id: str) -> Dict[str, Any]:
        """
        Make a request to the TED API.

        Args:
            endpoint: The API endpoint to call
            uniprot_id: The UniProt ID to query

        Returns:
            Dict[str, Any]: The response from the API
        """
        url = f"{self.TED_API_URL}{endpoint}/{uniprot_id}"
        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()

            if not response.content:
                return {}

            data = response.json()
            return data
        except requests.exceptions.HTTPError as e:
            logging.error(
                f"TED API request failed for {uniprot_id} with status code {response.status_code}: {e}"
            )
            raise
        except requests.exceptions.RequestException as e:
            logging.error(f"TED API request failed for {uniprot_id}: {e}")
            raise

    def get_TED_annotations(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Get TED annotations for a given UniProt ID.

        Args:
            uniprot_id: The UniProt ID to query

        Returns:
            List[Dict[str, Any]]: A list of TED annotations
        """
        endpoint = "summary"
        try:
            data = self._make_request(endpoint, uniprot_id)
        except (requests.exceptions.RequestException, ValueError) as e:
            logging.error(f"Error retrieving TED annotations for {uniprot_id}: {e}")
            return []

        if not data or "data" not in data:  # Check if the main 'data' key exists
            # Or simply rely on data.get('data', []) and subsequent checks
            logging.error(
                f"No 'data' field in API response for {uniprot_id} or response is empty."
            )
            return []

        annotations_raw = data.get("data", [])
        if not isinstance(annotations_raw, list):
            logging.error(f"Unexpected data format for {uniprot_id}: {annotations_raw}")
            return []

        processed_annotations = []

        for item in annotations_raw:
            if isinstance(item, dict):
                processed_annotations.append(
                    {
                        "uniprot_acc": uniprot_id,
                        "consensus_level": item.get("consensus_level"),
                        "chopping": item.get("chopping"),
                        "num_segments": item.get("num_segments"),
                        "cath_label": item.get("cath_label"),
                        "nres_domain": item.get("nres_domain"),
                    }
                )
        if not processed_annotations:
            logging.error(f"No valid annotations found for {uniprot_id}")
            return []
        logging.info(f"TED annotations for {uniprot_id}: {processed_annotations}")
        return processed_annotations


if __name__ == "__main__":
    client = TEDClient()

    # --- Test the TEDClient ---
    print("Testing TEDClient...")
    test_uniprot_id = "P12345"  # Replace with a valid UniProt ID for testing
    annotations = client.get_TED_annotations(test_uniprot_id)
    if annotations:
        print(f"TED annotations for {test_uniprot_id}:")
        for annotation in annotations:
            print(annotation)
    else:
        print(f"No TED annotations found for {test_uniprot_id}.")
    print("TEDClient test completed.")

    test_uniprot_invalid_id = (
        "INVALID_ID"  # Replace with an invalid UniProt ID for testing
    )
    annotations_invalid = client.get_TED_annotations(test_uniprot_invalid_id)
    if annotations_invalid:
        print(f"TED annotations for {test_uniprot_invalid_id}:")
        for annotation in annotations_invalid:
            print(annotation)
    else:
        print(f"No TED annotations found for {test_uniprot_invalid_id}.")
    print("TEDClient test with invalid ID completed.")
