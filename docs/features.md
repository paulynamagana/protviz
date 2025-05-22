---
layout: default
title: Features
nav_order: 2
description: "Discover the powerful features of Protviz for protein annotation visualisation."
---

# Features of Protviz
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
Protviz is equipped with a comprehensive suite of features designed to streamline the process of fetching, integrating, and visualising diverse protein annotations. Below is a more detailed look at what Protviz offers:
{: .fs-6 .fw-300 }

## 1. Extensive Data Retrieval Capabilities

Protviz connects to major bioinformatics databases, allowing you to gather a wide array of protein information with ease:

* **UniProt Integration**:
    * Quickly fetch the authoritative sequence length for any protein using its UniProt ID, forming the basis for all sequence-scaled visualisations.
* **PDBe (Protein Data Bank in Europe) Connectivity**:
    * **Structural Coverage**: Identify which regions of a protein have been experimentally determined and are available in PDB structures.
    * **Ligand Interactions**: Pinpoint specific residues involved in binding to various ligands, crucial for understanding protein function and drug design.
* **TED Database Access**:
    * Retrieve domain annotations from the TED (The Encyclopedia of Domains) database, often including CATH superfamily classifications.
* **AlphaFold Database (AFDB) Integration**:
    * **pLDDT Scores**: Fetch and visualise per-residue confidence scores (pLDDT) from AlphaFold predictions, offering insights into model quality.
    * **AlphaMissense Data**: Access AlphaMissense pathogenicity predictions to evaluate the potential impact of missense variants.
* **InterPro Database Support**:
    * Gain access to a wealth of domain and feature annotations from various InterPro member databases, such as:
        * **Pfam**: Visualise curated protein families and domains.
        * **CATH-Gene3D**: Display structural domain classifications.

## 2. Versatile Track-Based Visualisation

Protviz's core strength lies in its flexible track-based system, allowing you to layer multiple types of annotations onto a single, coherent plot:

* **`AxisTrack`**: The foundational track that displays the protein sequence as a horizontal axis, complete with numbered tick marks for easy reference. It automatically adapts to zoomed views.
* **`PDBTrack`**: Illustrates which parts of the protein are covered by PDB structures. Offers a `"full"` view (each PDB entry in a separate, labelled lane) or a `"collapse"` view (a merged representation of all covered regions).
* **`LigandInteractionTrack`**: Specifically designed to show residues that bind to ligands. Like `PDBTrack`, it supports `"full"` (per-ligand) and `"collapse"` (merged sites) plotting options.
* **`TEDDomainsTrack`**: Visualises domain segments from the TED database, with options to display each annotation uniquely or as a merged overview.
* **`AlphaFoldTrack`**: Dedicated to displaying AlphaFold data. It can simultaneously plot pLDDT scores (colour-coded by confidence) and average AlphaMissense pathogenicity predictions in distinct sub-tracks.
* **`InterProTrack`**: Displays domain signatures from InterPro member databases (e.g., Pfam, CATH-Gene3D). Supports `"full"` mode for detailed, per-signature rows with entry type and accession labels, and `"collapse"` mode for a merged view.
* **`CustomTrack`**: The most adaptable track, enabling you to plot your own data. You can define:
    * **Sequence Ranges**: Represent domains, motifs, or any region of interest as coloured bars.
    * **Point Features**: Mark specific residues (e.g., PTM sites, active site residues) using customisable markers.
    * Custom labels, colours, and grouping for each annotation.

## 3. Intuitive and Customisable Plotting

Protviz aims to make the creation of informative protein visualisations straightforward:

* **Unified Plotting Function**: The `plot_protein_tracks()` function serves as the central command to combine all your defined tracks into a final figure.
* **Zoom Functionality**: Easily focus on specific regions of interest by specifying `view_start_aa` and `view_end_aa` parameters, allowing for detailed examination of annotations in high-density areas.
* **Clear Labelling**: Automatic and customisable labelling for tracks, axes, and individual annotations ensures that your plots are easy to understand.
* **Save and Share**: Plots can be saved directly to image files (e.g., PNG) suitable for inclusion in presentations, publications, or web resources.
* **Aesthetic Control**: Options for colours, bar heights, spacing, and text sizes provide control over the visual appearance of the tracks and the overall plot.

By combining these features, Protviz provides a powerful yet accessible platform for creating rich, multi-layered visualisations of protein sequence annotations.