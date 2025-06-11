# src/protviz/data_retrieval/ted_client.py
import logging
import os
from typing import Any, Dict, List

import platformdirs
import requests
import requests_cache

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class TEDClient:
    """
    A class to interact with the TED annotations
    """

    TED_API_URL = "https://ted.cathdb.info/api/v1/uniprot/"

    def __init__(
        self,
        timeout: int = 10,
        cache_name: str = "ted_api_cache",
        expire_after: int = 86400,  # Cache expires after 1 day (86400 seconds)
        app_name: str = "protviz",
        app_author: str = "protviz_team",
    ):
        """
        Initialize the TED client with a timeout for requests.

        Args:
            timeout (int): Default timeout for API requests in seconds.
            cache_name (str): The base name for the cache file (e.g., 'pdbe_api_cache').
            expire_after (int): Time in seconds after which cached responses expire.
                                Use -1 for no expiration.
            app_name (str): The name of the application/package for platformdirs.
            app_author (str): The name for platformdirs.
        """
        self.timeout = timeout

        # --- Caching location
        # Get the platform-specific user cache directory
        user_cache_path = platformdirs.user_cache_dir(app_name, app_author)

        # Create the cache directory if it doesn't exist
        os.makedirs(user_cache_path, exist_ok=True)

        # construct the full cache file path
        cache_file_path = os.path.join(user_cache_path, f"{cache_name}.sqlite")
        logging.info(f"Using cache file at: {cache_file_path}")

        self.session = requests_cache.CachedSession(
            cache_file_path, backend="sqlite", expire_after=expire_after
        )

        self.session.headers.update({"Accept": "application/json"})

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

            if getattr(response, "from_cache", False):
                logging.info(f"Using cached response for {url}")
            else:
                logging.info(f"Making API request to {url}")

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

        except requests.exceptions.JSONDecodeError:
            cached = " (cached)" if getattr(response, "from_cache", False) else ""
            logging.error(
                f"TED API response for {url} is not valid JSON{cached}. Response content: {response.text[:200]}..."
            )
            raise ValueError(
                f"Invalid JSON response for {uniprot_id}. Check the API or the UniProt ID."
            )

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
            logging.info(f"No valid annotations found for {uniprot_id}")
            return []
        logging.info(f"TED annotations for {uniprot_id}: {processed_annotations}")
        return processed_annotations


if __name__ == "__main__":
    # Suggestion 3: Enhanced test block
    client = TEDClient()
    test_id_found = "Q53XC5"  # A protein kinase with known TED domains
    test_id_not_found = "A0A376RRE1"  # This one often doesn't have TED annotations

    print("\n--- RUN 1: POPULATING CACHE ---")

    # First call for an ID that should be found
    print(f"\n[1] Fetching annotations for {test_id_found}...")
    annotations_1 = client.get_TED_annotations(test_id_found)
    if annotations_1:
        print(f"Found {len(annotations_1)} annotations for {test_id_found}.")
        print(annotations_1[0])

    # First call for an ID that might not be found
    print(f"\n[1] Fetching annotations for {test_id_not_found}...")
    client.get_TED_annotations(test_id_not_found)

    print("\n\n--- RUN 2: HITTING CACHE ---")

    # Second call for the ID that was found
    print(f"\n[2] Fetching annotations for {test_id_found} again...")
    annotations_2 = client.get_TED_annotations(test_id_found)
    if annotations_2:
        print(
            f"Found {len(annotations_2)} annotations for {test_id_found} (from cache)."
        )

    # Second call for the ID that was not found
    print(f"\n[2] Fetching annotations for {test_id_not_found} again...")
    client.get_TED_annotations(test_id_not_found)

    print("\n--- Test complete ---")
