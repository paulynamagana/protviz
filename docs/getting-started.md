---
layout: default
title: Getting Started
nav_order: 3 
description: "A step-by-step guide to creating your first protein annotation plot with Protviz."
---

# Getting Started
{: .no_toc }

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

This guide will walk you through the fundamental steps to create your first protein annotation visualisation. We'll fetch some data for a chosen protein and display it using a few basic tracks.
{: .fs-6 .fw-300 }


# Basic Usage


The core workflow in Protviz involves these main stages:

1.  **Import necessary modules**: You'll need `plot_protein_tracks` for plotting, data retrieval clients (like `PDBeClient`, `AFDBClient`), and the specific track classes you wish to use (e.g., `AxisTrack`, `PDBTrack`).
2.  **Choose a UniProt ID**: This will be the protein you want to visualise.
3.  **Initialise data retrieval clients**: Create instances of the clients for the databases from which you want to fetch data.
4.  **Fetch protein information**:
    * Obtain the protein's sequence length using `get_protein_sequence_length()`. This is crucial for setting up the scale of your visualisation.
    * Employ the client instances to fetch specific annotation data (e.g., PDB coverage, AlphaFold data).
5.  **Create tracks**: For each piece of data you intend to visualise, instantiate a corresponding track object. Configure these tracks with the fetched data and any desired display options (labels, colours, plotting styles).
6.  **Plot the tracks**: Utilise the `plot_protein_tracks()` function, providing it with the protein ID, sequence length, and a list of the tracks you have created. You can also specify the figure width and an option to save the plot.


Here's a simple example of how to use Protviz to fetch data and plot tracks for a protein:
This example demonstrates a basic ProtViz workflow, fetching data from PDBe and UniProt and creating a simple plot.

```python

# 1. Import necessary modules
    from protviz import plot_protein_tracks
    from protviz.data_retrieval import PDBeClient, get_protein_sequence_length
    from protviz.tracks import AxisTrack, PDBTrack

    def main():
        # 2. Choose a UniProt ID
        uniprot_id = "P25494"  # Example: Human p53 protein
        # 3. Initialise Data Retrieval Clients
        pdbe_client = PDBeClient() 

        try:
            # 4. Fetch Protein Information
            seq_length = get_protein_sequence_length(uniprot_id) # Fetch sequence length
            print(f"Sequence length: {seq_length}")

            pdb_coverage = pdbe_client.get_pdb_coverage(uniprot_id) # Fetch PDB coverage data from PDBe
            print(f"PDB Coverage entries: {len(pdb_coverage)}")

            # 5. Create Tracks
            axis_track = AxisTrack(sequence_length=seq_length, label="Sequence") # The AxisTrack provides the main sequence axis
            pdb_track = PDBTrack(pdb_data=pdb_coverage, label="PDB Coverage") # The PDBTrack will show PDB coverage

            # 6. Plot the Tracks
            plot_protein_tracks(
                protein_id=uniprot_id,
                sequence_length=seq_length,
                tracks=[axis_track, pdb_track],
                figure_width=8, # Width of the output figure
                save_option=True # Saves the plot as <uniprot_id>_plot.png
            )

            print(f"Plot saved as {uniprot_id}_plot.png")

        except Exception as e:
            print(f"An error occurred: {e}")

        if __name__ == "__main__":
            main()
```

## Explanation of the code

* We commence by importing the necessary components from protviz.
* We define uniprot_id = "P25494" as our target protein.
* PDBeClient() is instantiated to interact with the PDBe database.
* get_protein_sequence_length(uniprot_id) fetches the total length of our protein.
* pdbe_client.get_pdb_coverage(uniprot_id) retrieves data concerning which parts of the protein are covered by PDB structures.
* We then create three types of tracks:
* AxisTrack: This is fundamental for displaying the sequence range.
    * PDBTrack: Configured to show a "collapsed" view of PDB coverage. This means overlapping regions from different PDB entries are merged into a single bar.
* Finally, plot_protein_tracks(...) takes all this information and generates the visualisation. The save_option=True argument ensures the plot is saved as O15245_plot.png.

After running this script, you will find an image file named O15245_plot.png in the same directory. This plot will display the sequence axis for O15245, a track indicating its PDB coverage.

This basic example demonstrates the power and simplicity of Protviz. In the following sections, we'll explore more advanced features, including different track types, data sources, and customisation options.