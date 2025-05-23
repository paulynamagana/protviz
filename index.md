---
title: Home
layout: home
nav_order: 1
description: ""
permalink: /
---
# **Protviz: Protein Annotation Visualizer**
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
![Language](https://img.shields.io/badge/Language-Python3-steelblue)
![Python Support](https://img.shields.io/badge/Python-3.9%20%7C3.10%20%7C%203.11%20%7C%203.12-blue)
![License](https://img.shields.io/badge/License-MIT-steelblue)
![Python Tests](https://github.com/paulynamagana/protviz/actions/workflows/python-tests.yml/badge.svg?branch=main)

ProtViz is a Python package designed to simplify the retrieval and visualization of protein annotations and structural information.
It provides a powerful and intuitive way to fetch data from multiple bioinformatics databases and plot this information along a protein sequence using a flexible track-based system.
{: .fs-6 .fw-300 }

## Key Features

Protviz simplifies complex protein data visualization with a range of powerful features:

* **Comprehensive Data Retrieval**:
    * Fetch protein sequence length directly from **UniProt**.
    * Retrieve PDB structure coverage and detailed ligand interaction data from the **Protein Data Bank in Europe (PDBe)**.
    * Get domain annotations from the **TED database**.
    * Fetch per-residue confidence scores (pLDDT) and AlphaMissense pathogenicity predictions from the **AlphaFold Database (AFDB)**.
    * Access domain information from **InterPro** member databases like Pfam and CATH-Gene3D.
<br>
<br>
* **Flexible Track-Based Visualization**: Protviz employs an intuitive track system to display diverse annotations:
    * **`AxisTrack`**: Essential for displaying the sequence axis and length context.
    * **`PDBTrack`**: Clearly shows regions covered by PDB structures, with options for detailed or collapsed views.
    * **`LigandInteractionTrack`**: Highlights residues involved in ligand binding.
    * **`TEDDomainsTrack`**: Visualizes domains as defined by the TED database.
    * **`AlphaFoldTrack`**: Displays AlphaFold's pLDDT and AlphaMissense scores along the sequence.
    * **`InterProTrack`**: Shows domain signatures from various InterPro member databases.
    * **`CustomTrack`**: Offers the ultimate flexibility to plot your own arbitrary annotations, whether they are sequence ranges or specific point features, with full control over appearance.
<br>
<br>
* **User-Friendly Plotting**:
    * Seamlessly combine multiple tracks into a single, coherent, and easy-to-interpret figure.
    * Effortlessly **zoom into specific regions** of the protein sequence for detailed analysis.
    * Conveniently save your plots in high quality for presentations and publications.

## Get Started Quickly

Jump right in by following our **[Getting Started](./docs/getting-started)** guide, or explore the **[Core Components](./docs/core-components/index)** to understand the data retrieval clients and visualization tracks in detail.

## Dependencies

Protviz relies on a few key Python libraries:

* numpy>=1.20
* matplotlib>=3.4
* requests>=2.20
* gemmi (specifically for parsing CIF files when fetching data from the AlphaFold Database)

We recommend installing Protviz in a dedicated virtual environment.

### License

ProtViz is distributed by an [MIT license](https://github.com/paulynamagana/protviz/tree/main/LICENSE)
