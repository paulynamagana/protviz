from protviz import plot_protein_tracks
from protviz.data_retrieval import PDBeClient, get_protein_sequence_length
from protviz.tracks import AxisTrack, LigandInteractionTrack, PDBTrack


def main():
    uniprot_id = (
        "P00533"  # EGFR - a receptor tyrosine kinase with many PDB entries and ligands
    )
    pdbe_client = PDBeClient()

    try:
        seq_length = get_protein_sequence_length(uniprot_id)
        print(f"Sequence length for {uniprot_id}: {seq_length}")

        # Fetch PDB coverage
        pdb_coverage_data = pdbe_client.get_pdb_coverage(uniprot_id)
        if pdb_coverage_data:
            print(f"Found {len(pdb_coverage_data)} PDB coverage entries.")
        else:
            print(f"No PDB coverage data found for {uniprot_id}.")

        # Fetch ligand interactions
        ligand_interaction_data = pdbe_client.get_pdb_ligand_interactions(uniprot_id)
        if ligand_interaction_data:
            print(f"Found {len(ligand_interaction_data)} ligand interaction contexts.")
        else:
            print(f"No ligand interaction data found for {uniprot_id}.")

        # Create tracks
        axis_trk = AxisTrack(sequence_length=seq_length, label="Sequence")

        pdb_trk_summary = PDBTrack(
            pdb_data=pdb_coverage_data,
            label="PDB Coverage",
            plotting_option="collapse",  # Show overall coverage
            color="skyblue",
        )

        ligand_trk_detailed = LigandInteractionTrack(
            interaction_data=ligand_interaction_data,
            label="Ligand Binding",
            plotting_option="full",  # Show each ligand separately
            show_ligand_labels=False,  # hide ligand labels
            site_height=0.15,  # Adjust height for clarity if many ligands
        )

        # Plot all tracks
        plot_protein_tracks(
            protein_id=uniprot_id,
            sequence_length=seq_length,
            tracks=[axis_trk, pdb_trk_summary, ligand_trk_detailed],
            figure_width=16,  # Adjust width as needed
            save_option=True,
        )
        print(f"PDBe example plot saved as {uniprot_id}_plot.png")

    except Exception as e:
        print(f"An error occurred during the PDBe example: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
