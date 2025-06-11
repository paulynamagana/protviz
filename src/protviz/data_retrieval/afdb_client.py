# src/protviz/data_retrieval/afdb_client.py
import csv
import io
import logging
import os
from typing import Any, Dict, List, Optional

import gemmi
import platformdirs
import requests
import requests_cache

# set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class AFDBClient:
    """
    A client to interact with the AlphaFold Database (AFDB) API
    and retrieve various annotation data.
    """

    AFDB_API_BASE_URL = "https://alphafold.ebi.ac.uk/api"

    def __init__(
        self,
        timeout: int = 20,
        cache_name: str = "afdb_api_cache",
        expire_after: int = 86400,  # Cache expires after 1 day (86400 seconds)
        app_name: str = "protviz",
        app_author: str = "protviz_team",
    ):
        """
        Initialize the AFDB client.
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

        # Construct the full cache file path
        cache_file_path = os.path.join(user_cache_path, f"{cache_name}.sqlite")
        logging.info(f"Using cache file at: {cache_file_path}")

        self.session = requests_cache.CachedSession(
            cache_file_path, backend="sqlite", expire_after=expire_after
        )

        self.session.headers.update({"Accept": "application/json"})

    def _make_api_request(self, endpoint_path: str) -> Any:
        """
        Makes a GET request to a specific AFDB API endpoint path.
        Expects JSON response.

        Args:
            endpoint_path (str): The API endpoint path (e.g., "/prediction/P12345").

        Returns:
            Any: The JSON response parsed into Python objects.

        Raises:
            requests.exceptions.HTTPError: If the request fails with a client/server error.
            requests.exceptions.RequestException: For other request-related issues.
            ValueError: If the response is not valid JSON.
        """
        url = f"{self.AFDB_API_BASE_URL}{endpoint_path}"
        logging.debug(f"Making AFDB API request to: {url}")
        try:
            response = self.session.get(
                url, timeout=self.timeout, headers={"Accept": "application/json"}
            )

            if getattr(response, "from_cache", False):
                logging.info(f"Using cached response for {url}")
            else:
                logging.info(f"Making API request to {url}")

            response.raise_for_status()

            if not response.content:
                logging.warning(f"Empty response content from {url}")
                return None

            return response.json()

        except requests.exceptions.HTTPError as e:
            logging.error(
                f"AFDB API HTTPError for {url}: {e} (Status: {e.response.status_code if e.response else 'N/A'})"
            )
            raise
        except requests.exceptions.JSONDecodeError:
            cached = " (cached)" if getattr(response, "from_cache", False) else ""
            logging.error(
                f"AFDB Failed to decode JSON from {url}{cached}. Response text: {response.text[:200]}..."
            )
            raise ValueError(
                f"Invalid JSON response from AFDB API at {url} at {endpoint_path}"
            )

        except requests.exceptions.RequestException as e:
            logging.error(f"AFDB API RequestException for {url}: {e}")
            raise

    def _fetch_file_content(self, file_url: str) -> str:
        """
        Fetches content from a given URL. Handles optional gzipping.

        Args:
            file_url (str): The URL of the file to fetch.

        Returns:
            str: The decoded content of the file as a string.

        Raises:
            requests.exceptions.RequestException: If the download fails.
            IOError: If decompression fails.
        """
        logging.debug(f"Fetching file from URL: {file_url}")
        try:
            response = self.session.get(file_url, timeout=self.timeout, stream=True)

            if getattr(response, "from_cache", False):
                logging.info(f"Using cached file content for {file_url}")
            else:
                logging.info(f"Fetching new file content from {file_url}")

            response.raise_for_status()

            content_bytes = response.content
            return content_bytes.decode("utf-8")
        except requests.exceptions.RequestException as e:
            logging.error(f"Failed to fetch file {file_url}: {e}")
            raise

    def get_alphafold_entry_details(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Gets the primary AlphaFold prediction entry details for a UniProt ID,
        which includes URLs to various data files.
        """
        endpoint = f"/prediction/{uniprot_id}"
        try:
            data = self._make_api_request(endpoint)
            if isinstance(data, list) and len(data) > 0:
                entry = data[0]
                logging.info(
                    f"Successfully fetched AlphaFold entry details for {uniprot_id}."
                )
                return {
                    "uniprotAccession": entry.get("uniprotAccession"),
                    "uniprotId": entry.get("uniprotId"),
                    "gene": entry.get("gene"),
                    "latestVersion": entry.get("latestVersion"),
                    "cifUrl": entry.get("cifUrl"),
                    "pdbUrl": entry.get("pdbUrl"),
                    "paeImageUrl": entry.get("paeImageUrl"),
                    "amAnnotationsUrl": entry.get("amAnnotationsUrl"),
                    "taxId": entry.get("taxId"),
                }
            else:
                logging.warning(
                    f"No AlphaFold prediction entries found for {uniprot_id} or unexpected format: {data}"
                )
                return None
        except (requests.exceptions.RequestException, ValueError):
            return None

    def _parse_plddt_from_cif_content(
        self, cif_content: str, uniprot_id: str
    ) -> List[Dict[str, Any]]:
        """
        Parses pLDDT scores from CIF file content string using gemmi.
        """
        logging.debug(f"Parsing pLDDT from CIF content for {uniprot_id}")
        plddt_scores = []
        try:
            # Read the CIF content from the string
            doc = gemmi.cif.read_string(cif_content)
            block = doc.sole_block()

            res_number_strs = block.find_values("_ma_qa_metric_local.label_seq_id")
            plddt_value_strs = block.find_values("_ma_qa_metric_local.metric_value")

            if not res_number_strs or not plddt_value_strs:
                logging.warning(
                    f"Could not find pLDDT tags or values in CIF for {uniprot_id}. "
                    f"Found {len(res_number_strs)} residue numbers and {len(plddt_value_strs)} pLDDT values."
                )
                return []

            if len(res_number_strs) != len(plddt_value_strs):
                logging.warning(
                    f"Mismatch in length of residue numbers ({len(res_number_strs)}) "
                    f"and pLDDT values ({len(plddt_value_strs)}) for {uniprot_id}. Cannot reliably parse."
                )
                return []

            for res_num_str, plddt_str in zip(res_number_strs, plddt_value_strs):
                try:
                    plddt_scores.append(
                        {"residue_number": int(res_num_str), "plddt": float(plddt_str)}
                    )
                except ValueError as e_val:
                    logging.warning(
                        f"Skipping pLDDT entry due to value conversion error for {uniprot_id}: "
                        f"res_num='{res_num_str}', plddt='{plddt_str}'. Error: {e_val}"
                    )
                    continue

            if not plddt_scores:
                logging.warning(
                    f"No pLDDT scores were successfully parsed from CIF for {uniprot_id}, "
                    "despite finding tags. Check data integrity or parsing logic for specific values."
                )

        except Exception as e:
            logging.error(
                f"Generic error parsing CIF content for pLDDT (UniProt: {uniprot_id}): {e}"
            )
            return []

        logging.info(
            f"Successfully parsed {len(plddt_scores)} pLDDT scores for {uniprot_id}"
        )
        return plddt_scores

    def _parse_alphamissense_csv_content(
        self, csv_content: str, uniprot_id: str
    ) -> List[Dict[str, Any]]:
        """Parses AlphaMissense scores from csv file content."""
        logging.debug(f"Parsing AlphaMissense csv content for {uniprot_id}")
        am_scores = []
        try:
            f = io.StringIO(csv_content)
            reader = csv.reader(f, delimiter=",")

            header = next(reader)  # Read the header row
            logging.debug(f"AlphaMissense header: {header}")

            expected_header_parts = ["protein_variant", "am_pathogenicity", "am_class"]
            if not all(part in header for part in expected_header_parts):
                logging.warning(
                    f"Unexpected AlphaMissense csv header for {uniprot_id}: {header}. Expected parts: {expected_header_parts}"
                )
                return []

            # Get column indices based on header
            try:
                variant_idx = header.index("protein_variant")
                pathogenicity_idx = header.index("am_pathogenicity")
                class_idx = header.index("am_class") if "am_class" in header else -1
            except ValueError as e_idx:
                logging.error(
                    f"Could not find required columns in AlphaMissense header for {uniprot_id}: {e_idx}"
                )
                return []

            for row_num, row in enumerate(reader):
                if not row or len(row) <= max(
                    variant_idx,
                    pathogenicity_idx,
                    (class_idx if class_idx != -1 else 0),
                ):
                    logging.warning(
                        f"Skipping malformed or short AlphaMissense row {row_num + 2} for {uniprot_id}: {row}"
                    )
                    continue

                variant_str = row[variant_idx]  # e.g., M1A
                try:
                    if (
                        not variant_str or len(variant_str) < 3
                    ):  # Basic validation for "M1A" format
                        logging.warning(
                            f"Invalid protein_variant format '{variant_str}' in row {row_num + 2} for {uniprot_id}"
                        )
                        continue

                    ref_aa = variant_str[0]
                    alt_aa = variant_str[-1]
                    res_num_str = variant_str[1:-1]

                    if not res_num_str.isdigit():
                        logging.warning(
                            f"Non-numeric residue number in variant '{variant_str}' in row {row_num + 2} for {uniprot_id}"
                        )
                        continue
                    res_num = int(res_num_str)

                    am_scores.append(
                        {
                            "residue_number": res_num,
                            "ref_aa": ref_aa,
                            "alt_aa": alt_aa,
                            "am_pathogenicity": float(row[pathogenicity_idx]),
                            "am_class": row[class_idx]
                            if class_idx != -1 and class_idx < len(row)
                            else "N/A",
                        }
                    )
                except (ValueError, IndexError, TypeError) as e_parse:
                    logging.warning(
                        f"Could not parse variant string or data in row {row_num + 2} ('{variant_str}') for {uniprot_id}: {e_parse}. Row: {row}"
                    )
                    continue
        except StopIteration:
            logging.warning(
                f"AlphaMissense file for {uniprot_id} seems empty or has no data rows after header."
            )
            return []
        except Exception as e:
            logging.error(f"Error parsing AlphaMissense csv for {uniprot_id}: {e}")
            return []

        logging.info(
            f"Successfully parsed {len(am_scores)} AlphaMissense scores for {uniprot_id}"
        )
        return am_scores

    def get_alphafold_data(
        self, uniprot_id: str, requested_data_types: List[str]
    ) -> Dict[str, Any]:
        """
        Fetches and processes various types of AlphaFold data for a UniProt ID.
        """
        results: Dict[str, Any] = {}
        logging.info(
            f"Fetching AlphaFold data for {uniprot_id}, requested types: {requested_data_types}"
        )

        entry_details = self.get_alphafold_entry_details(uniprot_id)
        if not entry_details:
            logging.warning(
                f"Could not get AlphaFold entry details for {uniprot_id}. Cannot fetch further data."
            )
            for req_type in requested_data_types:
                results[req_type] = []  # Ensure keys exist even if fetching fails
            return results

        if "plddt" in requested_data_types:
            cif_url = entry_details.get("cifUrl")
            if cif_url:
                try:
                    logging.info(f"Fetching CIF file for pLDDT: {cif_url}")
                    cif_content = self._fetch_file_content(cif_url)
                    if cif_content:
                        results["plddt"] = self._parse_plddt_from_cif_content(
                            cif_content, uniprot_id
                        )
                    else:
                        results["plddt"] = []
                except Exception as e:
                    logging.error(f"Failed to get pLDDT data for {uniprot_id}: {e}")
                    results["plddt"] = []
            else:
                logging.warning(
                    f"No cifUrl found for {uniprot_id}, cannot fetch pLDDT."
                )
                results["plddt"] = []

        if "alphamissense" in requested_data_types:
            am_url = entry_details.get("amAnnotationsUrl")
            tax_id = entry_details.get("taxId")
            if am_url:
                try:
                    logging.info(f"Fetching AlphaMissense file from: {am_url}")
                    am_content = self._fetch_file_content(am_url)
                    if am_content:
                        results["alphamissense"] = (
                            self._parse_alphamissense_csv_content(
                                am_content, uniprot_id
                            )
                        )
                    else:
                        results["alphamissense"] = []
                except Exception as e:
                    logging.error(
                        f"Failed to get AlphaMissense data for {uniprot_id}: {e}"
                    )
                    results["alphamissense"] = []
            else:
                if tax_id == 9606:
                    logging.warning(
                        f"No amAnnotationsUrl found for human entry {uniprot_id}."
                    )
                else:
                    logging.info(
                        f"No amAnnotationsUrl found for non-human entry {uniprot_id} (taxId: {tax_id})."
                    )
                results["alphamissense"] = []

        return results


if __name__ == "__main__":
    # Updated test block to demonstrate multi-level caching
    client = AFDBClient(timeout=30)

    # Use a human protein with AlphaMissense and a non-human one without
    human_id = "Q9BYF1"  # Nociceptin receptor, has AM
    non_human_id = "P0DTC2"  # Spike protein SARS-CoV-2

    # --- RUN 1: This will populate the cache for API calls and file downloads ---
    print("\n" + "=" * 40)
    print("      RUN 1: POPULATING CACHE")
    print("=" * 40)

    print(f"\n[1] Fetching pLDDT and AlphaMissense for HUMAN ({human_id})...")
    print("    (Expect 3 cache misses: API details, CIF download, CSV download)")
    data_human_1 = client.get_alphafold_data(
        human_id, requested_data_types=["plddt", "alphamissense"]
    )
    if data_human_1.get("plddt"):
        print(f"    -> Found {len(data_human_1['plddt'])} pLDDT scores.")
    if data_human_1.get("alphamissense"):
        print(
            f"    -> Found {len(data_human_1['alphamissense'])} AlphaMissense scores."
        )

    print(f"\n[1] Fetching pLDDT for NON-HUMAN ({non_human_id})...")
    print("    (Expect 2 cache misses: API details, CIF download)")
    data_non_human_1 = client.get_alphafold_data(
        non_human_id, requested_data_types=["plddt"]
    )
    if data_non_human_1.get("plddt"):
        print(f"    -> Found {len(data_non_human_1['plddt'])} pLDDT scores.")

    # --- RUN 2: This should hit the cache for all network requests ---
    print("\n\n" + "=" * 40)
    print("      RUN 2: HITTING CACHE")
    print("=" * 40)

    print(f"\n[2] Fetching pLDDT and AlphaMissense for HUMAN ({human_id}) again...")
    print("    (Expect 3 cache hits: API details, CIF download, CSV download)")
    data_human_2 = client.get_alphafold_data(
        human_id, requested_data_types=["plddt", "alphamissense"]
    )
    if data_human_2.get("plddt"):
        print(f"    -> Found {len(data_human_2['plddt'])} pLDDT scores (from cache).")
    if data_human_2.get("alphamissense"):
        print(
            f"    -> Found {len(data_human_2['alphamissense'])} AlphaMissense scores (from cache)."
        )

    print(f"\n[2] Fetching pLDDT for NON-HUMAN ({non_human_id}) again...")
    print("    (Expect 2 cache hits: API details, CIF download)")
    data_non_human_2 = client.get_alphafold_data(
        non_human_id, requested_data_types=["plddt"]
    )
    if data_non_human_2.get("plddt"):
        print(
            f"    -> Found {len(data_non_human_2['plddt'])} pLDDT scores (from cache)."
        )

    print("\nAFDBClient tests completed.")
