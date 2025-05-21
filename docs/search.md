---
title: Search
nav_order: 1
---

# Search
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

Just the Docs uses [lunr.js](https://lunrjs.com) to add a client-side search interface powered by a JSON index that Jekyll generates.
All search results are shown in an auto-complete style interface (there is no search results page).
By default, all generated HTML pages are indexed using the following data points:

- Page title
- Page content
- Page URL

## Enable search in configuration

In your site's `_config.yml`, enable search:

```python
Quick Start Code Example

This example demonstrates a basic ProtViz workflow, fetching data from PDBe and UniProt and creating a simple plot.

    from protviz import plot_protein_tracks
    from protviz.data_retrieval import PDBeClient, get_protein_sequence_length
    from protviz.tracks import AxisTrack, PDBTrack

    def main():
        uniprot_id = "P25494"  # Example: Human p53 protein
        pdbe_client = PDBeClient()

        try:
            seq_length = get_protein_sequence_length(uniprot_id)
            print(f"Sequence length: {seq_length}")

            pdb_coverage = pdbe_client.get_pdb_coverage(uniprot_id)
            print(f"PDB Coverage entries: {len(pdb_coverage)}")

            # Create tracks
            axis_track = AxisTrack(sequence_length=seq_length, label="Sequence")
            pdb_track = PDBTrack(pdb_data=pdb_coverage, label="PDB Coverage")

            # Plot tracks
            plot_protein_tracks(
                protein_id=uniprot_id,
                sequence_length=seq_length,
                tracks=[axis_track, pdb_track],
                figure_width=8
            )

        except Exception as e:
            print(f"An error occurred: {e}")

        if __name__ == "__main__":
            main()
```
