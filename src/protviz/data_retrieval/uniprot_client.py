# src/protviz/data_retrieval/uniprot_client.py
import requests


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
    headers = {"Acecept": "application/json"}

    try:
        response = requests.get(query, headers=headers)
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
    # Simple test
    try:
        # Example: Human Serum Albumin
        length_p02768 = get_protein_sequence_length("P02768")
        print(f"Length of P02768 (Human Serum Albumin): {length_p02768}")

        # Example: Spike protein of SARS-CoV-2
        length_p0dtc2 = get_protein_sequence_length("P0DTC2")
        print(f"Length of P0DTC2 (Spike Glycoprotein SARS-CoV-2): {length_p0dtc2}")

        # Example of an ID that might not exist to test error handling
        # get_protein_sequence_length("INVALIDID")
    except (requests.exceptions.HTTPError, KeyError, ValueError) as e:
        print(f"Error: {e}")
