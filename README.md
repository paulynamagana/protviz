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


## Installation

You need Python 3.9 or newer. To install Protviz, open a terminal and run:

```bash
pip install git+https://github.com/paulynamagana/protviz.git
```

## Quick start (no Python required)

Protviz comes with a `protviz` command, so you can make a figure without writing any code.
Give it a UniProt accession and it will fetch everything it can find and save a picture:

```bash
protviz P04637
```

That saves `P04637.png` in the folder you are in. A UniProt accession looks like `P04637`
or `Q9Y6K9` — it is not a gene name, so `TP53` will not work. You can find the accession
for your protein by searching its name at [uniprot.org](https://www.uniprot.org) and
copying the short code in the **Entry** column.

Some other things you can do:

```bash
protviz P04637 -o my_figure.png        # choose the file name
protviz P04637 --region 90-300         # zoom into amino acids 90 to 300
protviz P04637 --tracks pdb,pfam       # draw only the tracks you want
protviz P04637 --detail                # give every entry its own row
protviz P04637 --show                  # open a window instead of saving
protviz --list-tracks                  # see every track and what it shows
protviz --help                         # see all the options
```

If a database is unavailable or has nothing for your protein, Protviz says so and carries
on with the tracks it could fetch, rather than failing.

### Available tracks

| Track | Shows |
| --- | --- |
| `pdb` | Experimental structures from the PDB, as sequence coverage |
| `ligands` | Positions where a ligand binds the protein |
| `ted` | Structural domains predicted by TED |
| `pfam` | Protein families and domains from Pfam |
| `cath` | Structural domains from CATH-Gene3D |
| `plddt` | AlphaFold per-residue confidence |
| `alphamissense` | AlphaMissense average pathogenicity (slow — downloads a large file) |

`pdb`, `ligands`, `ted`, `pfam` and `plddt` are drawn by default. Use `--tracks all` to
include every one.

## Using Protviz from Python

The command line covers the common cases. Use the Python API when you want full control
over each track, or to add your own annotations with `CustomTrack`:

```python
from protviz import plot_protein_tracks
from protviz.data_retrieval import get_protein_sequence_length, PDBeClient, AFDBClient
from protviz.tracks import AxisTrack, PDBTrack, AlphaFoldTrack

uniprot_id = "O15245"

# Fetch the data you want
seq_length = get_protein_sequence_length(uniprot_id)
pdb_coverage = PDBeClient().get_pdb_coverage(uniprot_id)
alphafold_data = AFDBClient().get_alphafold_data(uniprot_id, requested_data_types=["plddt"])

# Build one track per kind of annotation
tracks = [
    AxisTrack(sequence_length=seq_length, label="Sequence"),
    PDBTrack(pdb_data=pdb_coverage, label="PDB", plotting_option="collapse"),
    AlphaFoldTrack(afdb_data=alphafold_data, plotting_options=["plddt"]),
]

# Draw them
plot_protein_tracks(
    protein_id=uniprot_id,
    sequence_length=seq_length,
    tracks=tracks,
    figure_width=12,
    save_path="my_figure.png",
    show=False,
)
```

`plot_protein_tracks` takes `save_path` to choose the output file and `show` to control
whether a window opens. Use `view_start_aa` and `view_end_aa` to zoom into a region.

## Dependencies

These are installed automatically with the package:

* numpy>=1.20
* matplotlib>=3.4
* requests>=2.32
* gemmi (for parsing CIF files from AlphaFold DB)
* requests-cache and platformdirs (for caching API responses between runs)

## Running Examples

The `examples/` folder contains scripts showing each data source in turn. Run one with:

```bash
python examples/example_afdb.py
```
