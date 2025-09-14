# Protviz: Protein Annotation Visualiser
![Python Support](https://img.shields.io/badge/Python-3.9%20%7C3.10%20%7C%203.11%20%7C%203.12-blue)
![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)
[![Python Tests](https://github.com/paulynamagana/protviz/actions/workflows/python-tests.yml/badge.svg?branch=main)](https://github.com/paulynamagana/protviz/actions/workflows/python-tests.yml)

Protviz is a Python package designed to retrieve and visualise various protein annotations and structural information. It allows users to fetch data from multiple bioinformatics databases and plot this information along a protein sequence using a flexible track-based system.

## Motivation

The goal of Protviz is to simplify the process of plotting protein annotations. This package is inspired by the [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html) library in R/Bioconductor, aiming to provide a similar, easy-to-use, track-based visualisation system for protein sequence data in Python.

I also wanted a way to plot data from resources but also be able to add custom annotations. Hope its helpful.



## Features

* **Data Retrieval**:
    * Fetch protein sequence length from UniProt.
    * Retrieve PDB coverage and ligand interaction data from PDBe.
    * Get TED domain annotations from the TED database.
    * Fetch pLDDT scores and AlphaMissense data from the AlphaFold Database (AFDB).

* **Track-Based Visualisation**:
    * **AxisTrack**: Displays the sequence axis with tick marks.
    * **PDBTrack**: Shows PDB structure coverage, with options to display as individual entries or a collapsed overview.
    * **LigandInteractionTrack**: Visualises ligand binding sites on the protein from PDB.
    * **TEDDomainsTrack**: Displays TED  annotations.
    * **AlphaFoldTrack**: Shows AlphaFold prediction metrics like pLDDT and average AlphaMissense pathogenicity scores.
    * **CustomTrack**: Allows plotting of arbitrary user-defined annotations (ranges or points) with customisable labels and colors.
    * **InterProTrack**: Displays InterPro annotations, like Pfam and CATH.


* **Core Plotting Functionality**:
    * Combines multiple tracks into a single, coherent plot.
    * Supports zooming into specific regions of the protein sequence.
    * Option to save plots to a file.


## Dependencies

The package requires the following Python libraries:

* numpy>=1.20
* matplotlib>=3.4
* requests>=2.20
* gemmi (for parsing CIF files from AlphaFold DB)

## Basic Usage

Here's a simple example of how to use Protviz to fetch data and plot tracks for a protein:

```python
from protviz import plot_protein_tracks
from protviz.data_retrieval import get_protein_sequence_length, PDBeClient, AFDBClient
from protviz.tracks import AxisTrack, PDBTrack, AlphaFoldTrack

def main():
    uniprot_id = "O15245"  # Example UniProt ID

    # Initialize clients
    pdbe_client = PDBeClient()
    afdb_client = AFDBClient()

    try:
        # Fetch sequence length
        seq_length = get_protein_sequence_length(uniprot_id)
        print(f"Sequence length for {uniprot_id}: {seq_length}")

        # Fetch PDB coverage data
        pdb_coverage = pdbe_client.get_pdb_coverage(uniprot_id)
        print(f"Found {len(pdb_coverage)} PDB coverage entries.")

        # Fetch AlphaFold data (pLDDT)
        alphafold_data = afdb_client.get_alphafold_data(
            uniprot_id,
            requested_data_types=["plddt"]
        )
        print(f"Fetched AlphaFold data. pLDDT entries: {len(alphafold_data.get('plddt', []))}")

        # Create tracks
        axis_track = AxisTrack(sequence_length=seq_length, label="Sequence")

        pdb_track_collapsed = PDBTrack(
            pdb_data=pdb_coverage,
            label="PDB (Collapsed)",
            plotting_option="collapse",
            bar_height=0.1,
            color="darkturquoise"
        )

        alphafold_track_plddt = AlphaFoldTrack(
            afdb_data=alphafold_data,
            plotting_options=["plddt"],
            main_label="AF pLDDT",
            sub_track_height=0.2
        )

        # Plot the tracks
        plot_protein_tracks(
            protein_id=uniprot_id,
            sequence_length=seq_length,
            tracks=[axis_track, pdb_track_collapsed, alphafold_track_plddt],
            figure_width=12,
            save_option=True # Saves the plot as <uniprot_id>_plot.png
        )

        print(f"Plot saved as {uniprot_id}_plot.png")

    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
```

## Running Examples

The package includes several example scripts (e.g., example_pdbe.py, example_afdb.py, example_ted.py, example_pdbe_zoom.py) in the root directory or an examples/ folder.

To run an example, navigate to the directory containing the scripts and execute it with Python:
```python
python examples/example_afdb.py
```
