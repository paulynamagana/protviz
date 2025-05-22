---
layout: default
title: Data Retrieval
parent: Core Components
nav_order: 1
description: "Understand the key building blocks of Protviz: Data Retrieval Clients and Visualisation Tracks."
---

# Core Components: Data Retrieval Clients and Visualisation Tracks
<button class="btn js-toggle-dark-mode">Dark mode</button>

<script>
const toggleDarkMode = document.querySelector('.js-toggle-dark-mode');

jtd.addEvent(toggleDarkMode, 'click', function(){
  if (jtd.getTheme() === 'dark') {
    jtd.setTheme('light');
    toggleDarkMode.textContent = 'Dark mode';
  } else {
    jtd.setTheme('dark');
    toggleDarkMode.textContent = 'Light mode';
  }
});
</script>
---
Protviz is built around two key types of components: **Data Retrieval Clients** that fetch information from various bioinformatics databases, and **Visualisation Tracks** that display this information. Understanding these components will help you leverage the full power of the package.
{: .fs-6 .fw-300 }

## Data Retrieval Clients

These clients are your interface to external bioinformatics databases, simplifying the process of fetching data for your protein of interest. Each client is tailored to a specific data source.

1.  `get_protein_sequence_length()`
    * **Purpose**: To quickly obtain the authoritative sequence length of a protein from UniProt. This is a fundamental piece of information required by many tracks to correctly scale their visualisations.
    * **Usage Example**:
        ```python
        from protviz.data_retrieval import get_protein_sequence_length

        uniprot_id = "P0DTC2"
        try:
            seq_len = get_protein_sequence_length(uniprot_id)
            print(f"The sequence length of {uniprot_id} is {seq_len}.")
        except ValueError as e:
            print(e)
        except Exception as e: # Other potential errors like network issues
            print(f"An unexpected error occurred: {e}")
        ```

2.  **`PDBeClient()`**
    * **Purpose**: To interact with the Protein Data Bank in Europe (PDBe) API, primarily for fetching information about PDB structures related to your UniProt ID, including structural coverage and ligand interactions.
    * **Initialisation and Usage Example**:
        ```python
        from protviz.data_retrieval import PDBeClient

        pdbe_client = PDBeClient()
        uniprot_id = "P00533" # EGFR, known for PDB entries and ligands

        pdb_structural_coverage = pdbe_client.get_pdb_coverage(uniprot_id)
        print(f"Found {len(pdb_structural_coverage)} PDB coverage entries for {uniprot_id}.")
        # Each item in pdb_structural_coverage is a dict like:
        # {'pdb_id': '2GS2', 'unp_start': 696, 'unp_end': 989}

        ligand_interaction_data = pdbe_client.get_pdb_ligand_interactions(uniprot_id)
        print(f"Found {len(ligand_interaction_data)} ligand interaction contexts for {uniprot_id}.")
        # Each item in ligand_interaction_data is a dict, e.g.:
        # {'ligand_id': 'STI', 'binding_site_uniprot_residues': [{'startIndex': 767, 'endIndex': 767, ...}, ...], 'pdb_id': '2ITP'}
        ```

3.  **`TEDClient()`**
    * **Purpose**: To retrieve domain annotations from the TED (The Encyclopedia of Domains) database. These annotations often include CATH domain classifications and segment information.
    * **Initialisation and Usage Example**:
        ```python
        from protviz.data_retrieval import TEDClient

        ted_client = TEDClient()
        uniprot_id = "O15245"

        ted_annotations = ted_client.get_TED_annotations(uniprot_id)
        print(f"Found {len(ted_annotations)} TED annotations for {uniprot_id}.")
        # Each item in ted_annotations is a dict like:
        # {'uniprot_acc': 'O15245', 'consensus_level': 'S100', 'chopping': '41-145_165-258', ...}
        ```

4.  **`AFDBClient()`**
    * **Purpose**: To connect with the AlphaFold Database (AFDB) and fetch prediction data, most notably pLDDT (per-residue confidence scores) and AlphaMissense pathogenicity predictions.
    * **Initialisation and Usage Example**:
        ```python
        from protviz.data_retrieval import AFDBClient

        afdb_client = AFDBClient()
        uniprot_id = "Q9BYF1" # Human protein with AlphaMissense data

        # You can request multiple data types
        alphafold_data_results = afdb_client.get_alphafold_data(
            uniprot_id,
            requested_data_types=["plddt", "alphamissense"]
        )

        plddt_scores = alphafold_data_results.get("plddt", [])
        print(f"Fetched {len(plddt_scores)} pLDDT scores for {uniprot_id}.")
        # Each item in plddt_scores is a dict like: {'residue_number': 1, 'plddt': 35.2}

        alphamissense_predictions = alphafold_data_results.get("alphamissense", [])
        print(f"Fetched {len(alphamissense_predictions)} AlphaMissense predictions for {uniprot_id}.")
        # Each item in alphamissense_predictions is a dict like:
        # {'residue_number': 1, 'ref_aa': 'M', 'alt_aa': 'A', 'am_pathogenicity': 0.02, 'am_class': 'likely_benign'}
        ```

5.  **`InterProClient()`**
    * **Purpose**: To fetch domain and feature annotations from various InterPro member databases (e.g., Pfam, CATH-Gene3D) by leveraging the general InterPro API.
    * **Initialisation and Usage Example**:
        ```python
        from protviz.data_retrieval import InterProClient

        interpro_client = InterProClient()
        uniprot_id = "P04637" # p53, has various InterPro entries

        pfam_domain_annotations = interpro_client.get_pfam_annotations(uniprot_id)
        print(f"Found {len(pfam_domain_annotations)} Pfam annotations for {uniprot_id}.")
        # Each item is a dict, e.g.:
        # {'accession': 'PF00870', 'name': 'P53 tumour suppressor family', 'description': 'P53 DNA-binding domain', ...}

        cath_gene3d_annotations = interpro_client.get_cathgene3d_annotations(uniprot_id)
        print(f"Found {len(cath_gene3d_annotations)} CATH-Gene3D annotations for {uniprot_id}.")
        # Similar structure to Pfam annotations, with CATH-Gene3D specific accessions.
        ```