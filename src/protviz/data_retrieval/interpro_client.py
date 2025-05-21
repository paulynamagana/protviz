# src/protviz/data_retrieval/interpro_client.py
import logging
from typing import Any, Dict, List, Optional

import requests

# Set up logging for this module
logger = logging.getLogger(__name__)


class InterProClient:
    """
    A client to interact with the InterPro API to retrieve protein domain annotations
    by parsing the general protein-centric InterPro entry list.
    """

    # This URL, when combined with "protein/uniprot/{uniprot_id}" as endpoint,
    # matches the user-provided API call structure.
    INTERPRO_API_BASE_URL = "https://www.ebi.ac.uk/interpro/api/entry/interpro"

    def __init__(self, timeout: int = 20):
        """
        Initialize the InterProClient.

        Args:
            timeout (int): Default timeout for API requests in seconds.
        """
        self.timeout = timeout
        self.session = requests.Session()

    def _fetch_protein_interpro_summary(
        self, uniprot_id: str
    ) -> Optional[Dict[str, Any]]:
        """
        Makes a GET request to fetch all InterPro entries associated with a UniProt ID.
        The typical endpoint for this is /api/protein/uniprot/{accession}, but the user's
        example uses /api/entry/interpro/protein/uniprot/{accession}. We will stick to user's pattern.

        Args:
            uniprot_id (str): The UniProt ID to query for.

        Returns:
            Optional[Dict[str, Any]]: The JSON response as a dictionary, or None if an error occurs.

        Raises:
            requests.exceptions.HTTPError: For 4xx/5xx errors if not handled as None.
            requests.exceptions.RequestException: For other network issues.
        """
        # This endpoint path, when appended to BASE_URL and then with uniprot_id,
        # forms: https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{uniprot_id}
        endpoint_path = "protein/uniprot"
        url = f"{self.INTERPRO_API_BASE_URL}/{endpoint_path}/{uniprot_id}"
        logger.debug(f"Making InterPro API request to: {url}")

        try:
            response = self.session.get(
                url, timeout=self.timeout, headers={"Accept": "application/json"}
            )
            response.raise_for_status()  # Raises an HTTPError for bad responses

            if not response.content:
                logger.warning(f"Empty response content from {url}")
                return None

            return response.json()

        except requests.exceptions.HTTPError as e:
            logger.error(
                f"InterPro API HTTPError for {url}: {e} (Status: {e.response.status_code if e.response else 'N/A'})"
            )
            if e.response is not None and e.response.status_code == 404:
                logger.info(
                    f"Resource not found (404) at {url} for UniProt ID {uniprot_id}"
                )
                return None
            raise  # Re-raise other HTTP errors
        except requests.exceptions.JSONDecodeError as e:
            logger.error(
                f"InterPro API JSONDecodeError for {url}: {e}. Response text: {response.text[:200]}"
            )
            return None
        except requests.exceptions.RequestException as e:
            logger.error(f"InterPro API RequestException for {url}: {e}")
            raise

    def _extract_member_db_annotations(
        self,
        uniprot_id: str,
        interpro_summary_data: Optional[Dict[str, Any]],
        member_db_key: str,
    ) -> List[Dict[str, Any]]:
        """
        Helper function to extract annotations for a specific member database (e.g., "pfam")
        from the general InterPro protein summary response.

        Args:
            uniprot_id (str): The UniProt ID queried.
            interpro_summary_data (Optional[Dict[str, Any]]): The parsed JSON response from
                                                             _fetch_protein_interpro_summary.
            member_db_key (str): The key for the member database (e.g., "pfam", "cathgene3d").

        Returns:
            List[Dict[str, Any]]: A list of processed annotations from the specified member database.
        """
        if not interpro_summary_data or not isinstance(
            interpro_summary_data.get("results"), list
        ):
            logger.warning(
                f"No 'results' found or invalid format in InterPro summary for {uniprot_id}."
            )
            return []

        annotations_list = []
        for interpro_entry in interpro_summary_data["results"]:
            if not isinstance(interpro_entry, dict):
                continue

            metadata = interpro_entry.get("metadata")
            if not isinstance(metadata, dict):
                continue

            parent_interpro_acc = metadata.get("accession")
            entry_type = metadata.get("type")
            name = metadata.get("name")
            member_databases = metadata.get("member_databases")

            if isinstance(member_databases, dict) and member_db_key in member_databases:
                db_signatures = member_databases[
                    member_db_key
                ]  # This is a dict of {db_id: description}
                if not isinstance(db_signatures, dict):
                    continue

                # Extract locations of the parent InterPro entry on this protein
                entry_locations_on_protein = []
                proteins_data = interpro_entry.get("proteins")
                if isinstance(proteins_data, list):
                    for protein_detail in proteins_data:
                        if (
                            isinstance(protein_detail, dict)
                            and protein_detail.get("accession", "").lower()
                            == uniprot_id.lower()
                            and isinstance(
                                protein_detail.get("entry_protein_locations"), list
                            )
                        ):
                            for loc_data in protein_detail["entry_protein_locations"]:
                                if isinstance(loc_data, dict) and isinstance(
                                    loc_data.get("fragments"), list
                                ):
                                    for fragment in loc_data["fragments"]:
                                        if (
                                            isinstance(fragment, dict)
                                            and fragment.get("start") is not None
                                            and fragment.get("end") is not None
                                        ):
                                            try:
                                                entry_locations_on_protein.append(
                                                    {
                                                        "start": int(fragment["start"]),
                                                        "end": int(fragment["end"]),
                                                    }
                                                )
                                            except ValueError:
                                                logger.warning(
                                                    f"Invalid start/end for fragment in {parent_interpro_acc} for {uniprot_id}"
                                                )

                if not entry_locations_on_protein:
                    logger.debug(
                        f"No locations found for InterPro entry {parent_interpro_acc} on protein {uniprot_id}."
                    )
                    # If InterPro entry has no locations, its member db signatures here also have no locations
                    # Or, you might decide to skip these signatures. For now, we'll include them if locations are found.
                    # continue # This would skip if parent IPR has no locations for the protein.

                for db_id, db_description in db_signatures.items():
                    annotations_list.append(
                        {
                            "accession": db_id,  # e.g., Pfam ID or CATH-Gene3D ID
                            "name": name,  # From the member_databases dict value
                            "description": db_description,
                            "entry_type": entry_type,
                            "interpro_accession": parent_interpro_acc,  # Parent InterPro entry
                            "locations": entry_locations_on_protein,  # Locations of the parent InterPro entry
                        }
                    )

        logger.info(
            f"Extracted {len(annotations_list)} {member_db_key} annotations for {uniprot_id} from InterPro summary."
        )
        logger.info(f"Annotations: {annotations_list}")
        return annotations_list

    def get_pfam_annotations(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Fetches Pfam domain annotations for a given UniProt ID by parsing the
        general InterPro protein summary.

        Args:
            uniprot_id (str): The UniProt ID to query.

        Returns:
            List[Dict[str, Any]]: A list of Pfam annotations.
        """
        try:
            interpro_summary_data = self._fetch_protein_interpro_summary(uniprot_id)
        except (requests.exceptions.RequestException, ValueError) as e:
            logger.error(
                f"Could not fetch InterPro summary for Pfam extraction ({uniprot_id}): {e}"
            )
            return []

        return self._extract_member_db_annotations(
            uniprot_id, interpro_summary_data, "pfam"
        )

    def get_cathgene3d_annotations(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Fetches CATH-Gene3D domain annotations for a given UniProt ID by parsing the
        general InterPro protein summary.

        Args:
            uniprot_id (str): The UniProt ID to query.

        Returns:
            List[Dict[str, Any]]: A list of CATH-Gene3D annotations.
        """
        try:
            interpro_summary_data = self._fetch_protein_interpro_summary(uniprot_id)
        except (requests.exceptions.RequestException, ValueError) as e:
            logger.error(
                f"Could not fetch InterPro summary for CATH-Gene3D extraction ({uniprot_id}): {e}"
            )
            return []

        return self._extract_member_db_annotations(
            uniprot_id, interpro_summary_data, "cathgene3d"
        )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(module)s - %(message)s",
    )

    client = InterProClient(timeout=25)

    test_uniprot_id = "P04637"  # p53 - Has Pfam and CATH-Gene3D via InterPro entries

    print(
        f"\n--- Testing Pfam Annotations for {test_uniprot_id} (from general InterPro summary) ---"
    )
    pfam_data = client.get_pfam_annotations(test_uniprot_id)
    if pfam_data:
        print(f"Found {len(pfam_data)} Pfam domain entries for {test_uniprot_id}:")
        for domain in pfam_data:
            print(
                f"  Pfam Acc: {domain['accession']}, Name: {domain['name']}, InterPro: {domain['interpro_accession']}, Type: {domain['entry_type']}"
            )
            print(f"    Locations: {domain['locations']}")
    else:
        print(f"No Pfam annotations found for {test_uniprot_id} or an error occurred.")

    print(
        f"\n--- Testing CATH-Gene3D Annotations for {test_uniprot_id} (from general InterPro summary) ---"
    )
    cath_data = client.get_cathgene3d_annotations(test_uniprot_id)
    if cath_data:
        print(
            f"Found {len(cath_data)} CATH-Gene3D domain entries for {test_uniprot_id}:"
        )
        for domain in cath_data:
            print(
                f"  CATH-Gene3D Acc: {domain['accession']}, Name: {domain['name']}, InterPro: {domain['interpro_accession']}, Type: {domain['entry_type']}"
            )
            print(f"    Locations: {domain['locations']}")
    else:
        print(
            f"No CATH-Gene3D annotations found for {test_uniprot_id} or an error occurred."
        )

    test_uniprot_id_no_data = "A0A000"  # Likely no data or 404
    print(
        f"\n--- Testing Pfam Annotations for {test_uniprot_id_no_data} (expecting none or error) ---"
    )
    pfam_data_none = client.get_pfam_annotations(test_uniprot_id_no_data)
    if not pfam_data_none:
        print(
            f"No Pfam annotations found for {test_uniprot_id_no_data} (as expected or error)."
        )
    else:
        print(
            f"Found {len(pfam_data_none)} Pfam domain entries for {test_uniprot_id_no_data} (unexpected)."
        )

    print("\nInterProClient tests completed.")
