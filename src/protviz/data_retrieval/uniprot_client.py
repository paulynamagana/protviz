# src/protviz/data_retrieval/uniprot_client.py
import logging
import os

import platformdirs
import requests
import requests_cache

# --- Caching Setup ---
app_name = "protviz"
app_author = "protviz_team"
cache_name = "uniprot_api_cache"
user_cache_path = platformdirs.user_cache_dir(app_name, app_author)
os.makedirs(user_cache_path, exist_ok=True)
cache_file_path = os.path.join(user_cache_path, f"{cache_name}.sqlite")

requests_cache.install_cache(
    cache_file_path, backend="sqlite", expire_after=86400
)  # 1 day expiration

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logging.info(f"Using UniProt cache file at: {cache_file_path}")


def get_protein_sequence_length(uniprot_id: str) -> int:
    """
    Fetch the protein sequence length from UniProt using the provided UniProt ID.

    Args:
    - uniprot_id (str): The UniProt ID of the protein.

    Returns:
    - int: The length of the protein sequence.

    Raises:
        requests.exceptions.HTTPError: If the request to UniProt fails.
        KeyError: If sequence information is not found in the response.
        ValueError: If the UniProt ID is invalid or protein not found.
    """
    UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb"
    query = f"{UNIPROT_API_URL}/{uniprot_id}"
    headers = {"Accept": "application/json"}

    try:
        response = requests.get(query, headers=headers)

        if getattr(response, "from_cache", False):
            logging.info(f"Using cached response for {query}")
        else:
            logging.info(f"Making API request to {query}")

        response.raise_for_status()  # Raise an error for bad responses

    except requests.exceptions.HTTPError:
        if response.status_code == 404:
            raise ValueError(f"UniProt ID '{uniprot_id}' not found.")
        else:
            raise

    data = response.json()

    if "sequence" in data and "length" in data["sequence"]:
        return int(data["sequence"]["length"])
    else:
        raise KeyError("Sequence information not found in the response.")


if __name__ == "__main__":
    # UPDATED: Test block to demonstrate caching

    # --- RUN 1: This will populate the cache ---
    print("\n" + "=" * 40)
    print("      RUN 1: POPULATING CACHE")
    print("=" * 40)

    try:
        print("\n[1] Fetching length for P02768 (Human Serum Albumin)...")
        length_p02768 = get_protein_sequence_length("P02768")
        print(f"    -> Length of P02768: {length_p02768}")

        print("\n[1] Fetching length for P0DTC2 (Spike Glycoprotein SARS-CoV-2)...")
        length_p0dtc2 = get_protein_sequence_length("P0DTC2")
        print(f"    -> Length of P0DTC2: {length_p0dtc2}")
    except (requests.exceptions.HTTPError, KeyError, ValueError) as e:
        print(f"Error during Run 1: {e}")

    # --- RUN 2: This should hit the cache ---
    print("\n\n" + "=" * 40)
    print("      RUN 2: HITTING CACHE")
    print("=" * 40)

    try:
        print("\n[2] Fetching length for P02768 again...")
        length_p02768_cached = get_protein_sequence_length("P02768")
        print(f"    -> Length of P02768: {length_p02768_cached} (from cache)")

        print("\n[2] Fetching length for an invalid ID to test error handling...")
        get_protein_sequence_length("INVALIDID")

    except (requests.exceptions.HTTPError, KeyError, ValueError) as e:
        print(f"    -> Successfully caught expected error: {e}")

    print("\nUniProt client tests completed.")
