from protviz import plot_protein_tracks
from protviz.data_retrieval import AFDBClient, PDBeClient, get_protein_sequence_length
from protviz.tracks import AlphaFoldTrack, AxisTrack, LigandInteractionTrack, PDBTrack


def main():
    uniprot_id = "O15245"
    pdbe_client = PDBeClient()
    afdb_client = AFDBClient()

    try:
        print(f"--- Protein: {uniprot_id} (for Ligand Track) ---")
        seq_length = get_protein_sequence_length(uniprot_id)
        print(f"Sequence length for {uniprot_id}: {seq_length}")

        print(f"\nFetching ligand interactions for {uniprot_id}...")
        ligand_data = pdbe_client.get_pdb_ligand_interactions(uniprot_id)
        if ligand_data:
            print(f"Found {len(ligand_data)} ligand interaction contexts.")
        else:
            print(f"No ligand interaction data found for {uniprot_id}.")

        # Fetch PDB coverage data also for this protein
        print(f"\nFetching PDB coverage for {uniprot_id}...")
        pdb_coverage = pdbe_client.get_pdb_coverage(uniprot_id)

        print(f"Fetching AlphaFold data for {uniprot_id} (pLDDT and AlphaMissense)...")
        # Request both data types from the client
        alphafold_data_fetched = afdb_client.get_alphafold_data(
            uniprot_id, requested_data_types=["plddt", "alphamissense"]
        )

        # Create tracks
        axis_trk = AxisTrack(
            sequence_length=seq_length,
            label="Sequence",
        )

        ligand_trk = LigandInteractionTrack(
            interaction_data=ligand_data,
            label="Ligands",
            show_ligand_labels=True,  # Show ligand IDs next to their interaction sites
            plotting_option="full",
        )

        pdb_trk_collapse = PDBTrack(
            pdb_data=pdb_coverage,
            label="PDB (Collapsed)",
            plotting_option="collapse",
            bar_height=0.1,
            color="darkturquoise",
        )

        # Create AlphaFoldTrack to plot both pLDDT and AlphaMissense
        alphafold_trk_both = AlphaFoldTrack(
            afdb_data=alphafold_data_fetched,
            plotting_options=["plddt", "alphamissense"],  # Specify what to plot
            main_label=" ",  # Optional main label for the track group
            plddt_label="pLDDT",  # Label for the pLDDT row
            alphamissense_label="AM",  # Label for the AlphaMissense row
            sub_track_height=0.05,  # Height of each row
            sub_track_spacing=0.03,  # Space between rows
        )

        #         # Create AlphaFoldTrack to plot only pLDDT
        #         alphafold_trk_plddt_only = AlphaFoldTrack(
        #             afdb_data=alphafold_data_fetched,
        #             plotting_options=["plddt"],
        #             main_label="AF pLDDT",
        #             sub_track_height=0.5
        #         )

        plot_protein_tracks(
            protein_id=uniprot_id,
            sequence_length=seq_length,
            tracks=[axis_trk, pdb_trk_collapse, ligand_trk, alphafold_trk_both],
            figure_width=16,
        )

        print("Done with ligand example.")

    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
