Basic Usage
===========

Here's a simple example of how to use Protviz to fetch data and plot tracks for a protein::

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

Running Examples
----------------

The package includes several example scripts (e.g., example_pdbe.py, example_afdb.py, example_ted.py, example_pdbe_zoom.py) in the root directory or an examples/ folder.

To run an example, navigate to the directory containing the scripts and execute it with Python::

    python examples/example_afdb.py

Dependencies
------------

The package requires the following Python libraries:

* numpy>=1.20
* matplotlib>=3.4
* requests>=2.20
* gemmi (for parsing CIF files from AlphaFold DB)