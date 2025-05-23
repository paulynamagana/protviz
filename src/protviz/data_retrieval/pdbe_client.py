# src/protviz/data_retrieval/pdbe_client.py
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


class PDBeClient:
    """
    Client for fetching data from the Protein Data Bank in Europe (PDBe) APIs,
    with built-in, user-directory-based caching.
    """

    PDBE_GRAPH_API_BASE_URL = "https://www.ebi.ac.uk/pdbe/graph-api"

    def __init__(
        self,
        timeout: int = 10,
        cache_name: str = "pdbe_api_cache",
        expire_after: int = 86400,  # Cache expires after 1 day (86400 seconds)
        app_name: str = "protviz",  # Your package/app name
        app_author: str = "ProtvizTeam",  # Your team/author name
    ):
        """
        Initialises the PDBeClient with caching enabled in the user cache dir.

        Args:
            timeout (int): Default timeout for API requests in seconds.
            cache_name (str): The base name for the cache file (e.g., 'pdbe_api_cache').
            expire_after (int): Time in seconds after which cached responses expire.
                                Use -1 for no expiration.
            app_name (str): The name of your application/package for platformdirs.
            app_author (str): The author/vendor name for platformdirs.
        """
        self.timeout = timeout

        # --- Caching Location Improved Here ---
        # Get the platform-specific user cache directory
        user_cache_path = platformdirs.user_cache_dir(app_name, app_author)

        # Create the directory if it doesn't exist
        os.makedirs(user_cache_path, exist_ok=True)

        # Construct the full path to the cache database file
        cache_file_path = os.path.join(user_cache_path, f"{cache_name}.sqlite")
        logging.info(f"Using cache file at: {cache_file_path}")

        # Create a CachedSession using the full path
        self.session = requests_cache.CachedSession(
            cache_file_path,  # Use the full path now
            backend="sqlite",
            expire_after=expire_after,
        )
        # --- End Caching Location Improvement ---

        self.session.headers.update({"Accept": "application/json"})

    def _make_request(self, endpoint: str, uniprot_id: str) -> Dict[str, Any]:
        """
        Helper method to make a GET request to a PDBe graph API endpoint,
        utilizing the cached session.
        """
        url = f"{self.PDBE_GRAPH_API_BASE_URL}{endpoint}/{uniprot_id}"
        try:
            response = self.session.get(url, timeout=self.timeout)

            if getattr(response, "from_cache", False):
                logging.info(f"Cache hit for {url}")
            else:
                logging.info(f"Cache miss for {url}. Fetched from API.")

            response.raise_for_status()

            if not response.content:
                return {}

            data = response.json()
            return data

        except requests.exceptions.HTTPError as e:
            logging.error(
                f"PDBe API request failed for {url} with status {response.status_code}: {e}"
            )
            raise
        except requests.exceptions.JSONDecodeError:
            cached = " (cached)" if getattr(response, "from_cache", False) else ""
            logging.error(
                f"Failed to decode JSON from {url}{cached}. Response text: {response.text[:200]}..."
            )
            raise ValueError(
                f"Invalid JSON response from PDBe API for {uniprot_id} at {endpoint}."
            )
        except requests.exceptions.RequestException as e:
            logging.error(f"PDBe API request failed for {url}: {e}")
            raise

    def get_pdb_coverage(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Fetch the PDB structure coverage for a given UniProt ID from PDBe.
        Uses the "best_structures" mapping endpoint.

        Args:
            uniprot_id (str): The UniProt ID of the protein.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries, each containing PDB coverage information
                                  (pdb_id, unp_start, unp_end).
                                  Returns an empty list if no coverage is found.
        """
        endpoint = "/mappings/best_structures"
        try:
            data = self._make_request(endpoint, uniprot_id)
        except (requests.exceptions.RequestException, ValueError) as e:
            logging.error(f"Could not fetch PDB coverage for {uniprot_id}: {e}")
            return []  # Return empty on error to allow downstream processing

        if not data or uniprot_id not in data:
            return []

        pdb_entries_raw = data.get(uniprot_id, [])
        processed_entries = []
        for entry in pdb_entries_raw:
            if "pdb_id" in entry and "unp_start" in entry and "unp_end" in entry:
                processed_entries.append(
                    {
                        "pdb_id": entry["pdb_id"],
                        "unp_start": int(entry["unp_start"]),
                        "unp_end": int(entry["unp_end"]),
                    }
                )
            else:
                logging.warning(
                    f"Warning: Skipping PDB coverage entry due to missing fields for {uniprot_id}: {entry}"
                )
        logging.info(f"All PDB entries found for {uniprot_id}: {processed_entries}")
        return processed_entries

    def get_pdb_ligand_interactions(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Fetch residues interacting with ligands for a given UniProt ID from PDBe.

        Args:
            uniprot_id (str): The UniProt ID of the protein.

        Returns:
            List[Dict[str, Any]]: A list of dictionaries, where each dictionary
                                  represents a ligand and its interacting residues.
                                  Includes:
                                  - 'ligand_accession': Accession ID of the ligand (e.g., "VIB")
                                  - 'interacting_residues': List of dicts, each with:
                                    - 'startIndex': UniProt start residue index of interaction
                                    - 'endIndex': UniProt end residue index of interaction
                                    - 'startCode': Amino acid code at startIndex
                                    - 'endCode': Amino acid code at endIndex
                                  Returns an empty list if no relevant data is found or on error.
        """
        endpoint = "/uniprot/ligand_sites"  # endpoint for ligand interactions
        try:
            data = self._make_request(endpoint, uniprot_id)

        except (requests.exceptions.RequestException, ValueError) as e:
            logging.error(f"Could not fetch ligand interactions for {uniprot_id}: {e}")
            return []

        if not data or uniprot_id not in data:
            return []

        protein_data = data.get(uniprot_id)
        if not isinstance(protein_data, dict):
            logging.error(
                f"Unexpected data format for ligand interactions for {uniprot_id}: {protein_data}"
            )
            return []

        if protein_data.get("dataType") != "LIGAND BINDING SITES":
            logging.error(
                f"Unexpected data type for ligand interactions for {uniprot_id}: {protein_data.get('dataType')}"
            )
            return []

        ligand_interactions_raw = protein_data.get("data", [])
        if not isinstance(ligand_interactions_raw, list):
            logging.error(f"Ligands data for {uniprot_id} is not a list")
            return []

        processed_interactions = []

        for (
            ligand_entry
        ) in ligand_interactions_raw:  # Iterating over list of PDB entries with ligands
            if not isinstance(ligand_entry, dict):
                logging.error(
                    f"Unexpected ligand entry format for {uniprot_id}: {ligand_entry}"
                )
                continue

            ligand_accession = ligand_entry.get("accession")  # Or 'chemical_name'

            if not ligand_accession:
                logging.warning(
                    f"Warning: Skipping ligand entry due to missing accession for {uniprot_id}: {ligand_entry}"
                )
                continue

            residues_raw = ligand_entry.get("residues", [])
            if not isinstance(residues_raw, list):
                logging.error(
                    f"Unexpected residues format for {uniprot_id}: {residues_raw}"
                )
                continue

            interacting_residues_details = []
            for res_info in residues_raw:
                if not isinstance(res_info, dict):
                    logging.error(
                        f"Skipping invalid residue info for ligand {ligand_accession}."
                    )
                    continue

                start_index = res_info.get("startIndex")
                end_index = res_info.get("endIndex")
                start_code = res_info.get("startCode")
                end_code = res_info.get("endCode")

                if start_index is not None and end_index is not None:
                    interacting_residues_details.append(
                        {
                            "startIndex": int(start_index),
                            "endIndex": int(end_index),
                            "startCode": start_code,
                            "endCode": end_code,
                        }
                    )

            if interacting_residues_details:
                processed_interactions.append(
                    {
                        "ligand_id": ligand_accession,  # Unique identifier for the ligand
                        "binding_site_uniprot_residues": interacting_residues_details,
                        "pdb_id": ligand_entry.get("pdb_id", ""),
                    }
                )
        logging.info(
            f"All ligand interactions found for {uniprot_id}: {processed_interactions}"
        )
        return processed_interactions

    def clear_cache(self):
        """Clears the entire requests cache."""
        logging.info("Clearing PDBe API cache...")
        self.session.cache.clear()
        logging.info("Cache cleared.")


if __name__ == "__main__":
    client = PDBeClient()  # Create an instance of the client

    # Optionally clear cache before starting a fresh test run
    # client.clear_cache()

    # Define the IDs to test
    uniprot_id_with_coverage = "P07550"  # TGF-beta receptor type-1
    uniprot_id_without_coverage = "A0A023GPI8"  # Might have no "best_structures"
    uniprot_id_with_ligands = "P00533"  # EGFR
    uniprot_id_few_ligands = "P07550"  # Same as with_coverage
    uniprot_id_no_ligands_expected = "A0A023GPI8"  # Same as without_coverage

    print("=" * 30)
    print("      FIRST RUN (Populating Cache)     ")
    print("=" * 30)

    # --- Test PDB Coverage (Run 1) ---
    print("\n--- Testing PDB Coverage (Run 1) ---")
    print(f"Fetching PDB coverage for {uniprot_id_with_coverage}...")
    client.get_pdb_coverage(uniprot_id_with_coverage)
    print(f"Fetching PDB coverage for {uniprot_id_without_coverage}...")
    client.get_pdb_coverage(uniprot_id_without_coverage)

    # --- Test Ligand Interactions (Run 1) ---
    print("\n--- Testing Ligand Interactions (Run 1) ---")
    print(f"Fetching ligand interactions for {uniprot_id_with_ligands}...")
    client.get_pdb_ligand_interactions(uniprot_id_with_ligands)
    print(f"Fetching ligand interactions for {uniprot_id_few_ligands}...")
    client.get_pdb_ligand_interactions(uniprot_id_few_ligands)
    print(f"Fetching ligand interactions for {uniprot_id_no_ligands_expected}...")
    client.get_pdb_ligand_interactions(uniprot_id_no_ligands_expected)

    print("\n\n")
    print("=" * 30)
    print("      SECOND RUN (Testing Cache)       ")
    print("=" * 30)

    # --- Test PDB Coverage (Run 2 - Should Hit Cache) ---
    print("\n--- Testing PDB Coverage (Run 2) ---")
    print(f"Fetching PDB coverage for {uniprot_id_with_coverage}...")
    coverage_data = client.get_pdb_coverage(uniprot_id_with_coverage)
    if coverage_data:
        print(f"Found {len(coverage_data)} PDB entries.")
    else:
        print("No PDB coverage found.")

    print(f"\nFetching PDB coverage for {uniprot_id_without_coverage}...")
    coverage_data_none = client.get_pdb_coverage(uniprot_id_without_coverage)
    if coverage_data_none:
        print(f"Found {len(coverage_data_none)} PDB entries.")
    else:
        print("No PDB coverage found.")

    # --- Test Ligand Interactions (Run 2 - Should Hit Cache) ---
    print("\n--- Testing Ligand Interactions (Run 2) ---")
    print(f"Fetching ligand interactions for {uniprot_id_with_ligands}...")
    ligand_data = client.get_pdb_ligand_interactions(uniprot_id_with_ligands)
    if ligand_data:
        print(f"Found {len(ligand_data)} ligand contexts.")
    else:
        print("No ligand interactions found.")

    print(f"\nFetching ligand interactions for {uniprot_id_few_ligands}...")
    ligand_data_tgf = client.get_pdb_ligand_interactions(uniprot_id_few_ligands)
    if ligand_data_tgf:
        print(f"Found {len(ligand_data_tgf)} ligand contexts.")
    else:
        print("No ligand interactions found.")

    print(f"\nFetching ligand interactions for {uniprot_id_no_ligands_expected}...")
    ligand_data_none = client.get_pdb_ligand_interactions(
        uniprot_id_no_ligands_expected
    )
    if not ligand_data_none:
        print("No ligand interactions found.")
    else:
        print(f"Found {len(ligand_data_none)} ligand contexts.")

    print("\n--- Cache Test Complete ---")
