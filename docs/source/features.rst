Features
========

Protviz is a Python package designed to retrieve and visualize various protein annotations and structural information. It allows users to fetch data from multiple bioinformatics databases and plot this information along a protein sequence using a flexible track-based system.

Data Retrieval
--------------

* Fetch protein sequence length from UniProt.
* Retrieve PDB coverage and ligand interaction data from PDBe.
* Get TED domain annotations from the TED database.
* Fetch pLDDT scores and AlphaMissense data from the AlphaFold Database (AFDB).

Track-Based Visualization
-------------------------

* **AxisTrack**: Displays the sequence axis with tick marks.
* **PDBTrack**: Shows PDB structure coverage, with options to display as individual entries or a collapsed overview.
* **LigandInteractionTrack**: Visualizes ligand binding sites on the protein.
* **TEDDomainsTrack**: Displays TED domain annotations.
* **AlphaFoldTrack**: Shows AlphaFold prediction metrics like pLDDT and average AlphaMissense pathogenicity scores.
* **CustomTrack**: Allows plotting of arbitrary user-defined annotations (ranges or points) with customizable labels and colors.

Core Plotting Functionality
--------------------------

* Combines multiple tracks into a single, coherent plot.
* Supports zooming into specific regions of the protein sequence.
* Option to save plots to a file.