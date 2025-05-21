# src/protviz/data_retrieval/pdbe_client.py
import logging
from typing import Any, Dict, List

import requests

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class PDBeClient:
    """
    Client for fetching data from the Protein Data Bank in Europe (PDBe) APIs.
    """

    PDBE_GRAPH_API_BASE_URL = "https://www.ebi.ac.uk/pdbe/graph-api"

    def __init__(self, timeout: int = 10):
        """
        Initialises the PDBeClient.

        Args:
            timeout (int): Default timeout for API requests in seconds.
        """
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})

    def _make_request(self, endpoint: str, uniprot_id: str) -> Dict[str, Any]:
        """
        Helper method to make a GET request to a PDBe graph API endpoint.

        Args:
            endpoint (str): The specific API endpoint path (e.g., "/mappings/best_structures").
            uniprot_id (str): The UniProt ID to query for.

        Returns:
            Dict[str, Any]: The JSON response as a dictionary.

        Raises:
            requests.exceptions.HTTPError: If the request to PDBe fails with a client/server error.
            requests.exceptions.RequestException: For other request-related issues (timeout, connection error).
            ValueError: If the response is not valid JSON or the UniProt ID is not in the response.
        """
        url = f"{self.PDBE_GRAPH_API_BASE_URL}{endpoint}/{uniprot_id}"
        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()  # Raises an HTTPError for bad responses (4XX or 5XX)

            # Handle cases where API returns 200 OK with an empty body or non-JSON content
            if not response.content:
                # If UniProt ID not found, some graph endpoints return 200 OK with {"uniprot_id": []}
                # or simply an empty dict {} for truly non-existent/malformed top-level queries.
                # This check handles if the response is completely empty.
                return {}

            data = response.json()
            return data

        except requests.exceptions.HTTPError as e:
            # Specific handling for 404 could be nuanced depending on the endpoint
            # For graph API, a 404 might mean the UniProt ID format is wrong or the base endpoint path is wrong.
            # If the UniProt ID is valid but has no data, it often returns 200 OK with an empty list for that ID.
            logging.error(
                f"PDBe API request failed for {url} with status {response.status_code}: {e}"
            )
            raise
        except requests.exceptions.JSONDecodeError:
            logging.error(
                f"Failed to decode JSON from {url}. Response text: {response.text[:200]}..."
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


if __name__ == "__main__":
    client = PDBeClient()  # Create an instance of the client

    # --- Test get_pdb_coverage ---
    print("--- Testing PDB Coverage ---")
    uniprot_id_with_coverage = "P07550"  # TGF-beta receptor type-1
    print(f"\nFetching PDB coverage for {uniprot_id_with_coverage}...")
    coverage_data = client.get_pdb_coverage(uniprot_id_with_coverage)
    if coverage_data:
        print(f"Found {len(coverage_data)} PDB entries for {uniprot_id_with_coverage}:")
        for entry in coverage_data[:3]:  # Print first 3
            print(
                f"  PDB: {entry['pdb_id']}, UniProt Range: {entry['unp_start']}-{entry['unp_end']}"
            )
    else:
        print(f"No PDB coverage found for {uniprot_id_with_coverage}.")

    uniprot_id_without_coverage = "A0A023GPI8"  # Might have no "best_structures"
    print(f"\nFetching PDB coverage for {uniprot_id_without_coverage}...")
    coverage_data_none = client.get_pdb_coverage(uniprot_id_without_coverage)
    if coverage_data_none:
        print(
            f"Found {len(coverage_data_none)} PDB entries for {uniprot_id_without_coverage}."
        )
    else:
        print(f"No PDB coverage found for {uniprot_id_without_coverage}.")

    # --- Test get_pdb_ligand_interactions ---
    print("\n--- Testing Ligand Interactions ---")
    uniprot_id_with_ligands = "P00533"  # EGFR, known to have many PDBs with ligands
    print(f"\nFetching ligand interactions for {uniprot_id_with_ligands}...")
    ligand_data = client.get_pdb_ligand_interactions(uniprot_id_with_ligands)
    if ligand_data:
        print(
            f"Found {len(ligand_data)} ligand interaction contexts for {uniprot_id_with_ligands}:"
        )
        for interaction in ligand_data[:3]:  # Print first 3 contexts
            print(f"  PDB: {interaction['pdb_id']}, Ligand: {interaction['ligand_id']}")
            print(
                f"    Interacting UniProt Residues: {interaction['binding_site_uniprot_residues'][:5]}..."
            )  # First 5 residues
            # print(f"    Detailed Residues: {interaction['interacting_residues_details'][:2]}") # First 2 detailed
    else:
        print(f"No ligand interactions found for {uniprot_id_with_ligands}.")

    uniprot_id_few_ligands = (
        "P07550"  # Check TGF-beta receptor type-1 again for ligands
    )
    print(f"\nFetching ligand interactions for {uniprot_id_few_ligands}...")
    ligand_data_tgf = client.get_pdb_ligand_interactions(uniprot_id_few_ligands)
    if ligand_data_tgf:
        print(
            f"Found {len(ligand_data_tgf)} ligand interaction contexts for {uniprot_id_few_ligands}:"
        )
        for interaction in ligand_data_tgf:
            print(f"  PDB: {interaction['pdb_id']}, Ligand: {interaction['ligand_id']}")
            print(
                f"    Interacting UniProt Residues: {interaction['binding_site_uniprot_residues']}"
            )
    else:
        print(f"No ligand interactions found for {uniprot_id_few_ligands}.")

    # Test with an ID that might not have ligand data from this specific endpoint
    uniprot_id_no_ligands_expected = "A0A023GPI8"
    print(f"\nFetching ligand interactions for {uniprot_id_no_ligands_expected}...")
    ligand_data_none = client.get_pdb_ligand_interactions(
        uniprot_id_no_ligands_expected
    )
    if not ligand_data_none:
        print(
            f"No ligand interactions found for {uniprot_id_no_ligands_expected} (as expected or API has no data)."
        )
