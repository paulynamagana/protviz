from protviz import plot_protein_tracks
from protviz.data_retrieval import AFDBClient, PDBeClient, get_protein_sequence_length
from protviz.tracks import (
    AlphaFoldTrack,
    AxisTrack,
    CustomTrack,
    LigandInteractionTrack,
    PDBTrack,
)


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

        print(f"\nFetching PDB coverage for {uniprot_id}...")
        pdbe_client = PDBeClient()  # Create an instance of the PDBe client
        pdb_coverage_data = pdbe_client.get_pdb_coverage(uniprot_id)
        if pdb_coverage_data:
            print(f"Found {len(pdb_coverage_data)} PDB entries.")
        else:
            print("No PDB coverage data found.")

        print(f"Fetching AlphaFold data for {uniprot_id} (pLDDT and AlphaMissense)...")
        # Request both data types from the client
        alphafold_data_fetched = afdb_client.get_alphafold_data(
            uniprot_id, requested_data_types=["plddt", "alphamissense"]
        )

        # --- Zoomed View ---
        zoom_s, zoom_e = 150, 300
        if seq_length < zoom_e:  # Adjust zoom if protein is shorter
            zoom_s = seq_length // 4
            zoom_e = (seq_length // 4) * 3
        if zoom_s >= zoom_e:
            zoom_s = 1
            zoom_e = seq_length

        print(f"\nPlotting Zoomed View ({zoom_s}-{zoom_e})...")
        axis_trk_zoom = AxisTrack(
            sequence_length=seq_length,
            label="Sequence",  # Pass total length
        )

        ptm_annotations = [
            {
                "position": 245,
                "label": "S45-P",
                "color": "blue",
                "display_type": "marker",
                "marker_symbol": "P",
            },
            {
                "position": 188,
                "label": "Y88-P",
                "color": "#FF00FF",
                "display_type": "marker",
                "marker_symbol": "*",
            },
            {"position": 150, "label": "T150-Gly", "color": "green"},
            {
                "position": 200,
                "label": "Important",
                "color": "red",
                "display_type": "marker",
                "marker_symbol": "*",
            },
        ]
        ptms_track = CustomTrack(
            annotation_data=ptm_annotations,
            label="Modifications",  # Main track label
            ann_height=0.05,  # Height for the single row of PTMs
            show_row_labels=True,
            show_ann_labels=True,
            padding=0.05,
        )

        pdb_trk_detail_zoom = PDBTrack(
            pdb_data=pdb_coverage_data,
            label="PDB",
            plotting_option="collapse",
            color="lightcoral",
        )

        ligand_trk = LigandInteractionTrack(
            interaction_data=ligand_data,
            label="Ligand \n Binding Sites",
            show_ligand_labels=True,
            plotting_option="collapse",
            site_height=0.1,
        )

        # Create AlphaFoldTrack to plot both pLDDT and AlphaMissense
        alphafold_trk_both = AlphaFoldTrack(
            afdb_data=alphafold_data_fetched,
            plotting_options=["plddt"],  # Specify what to plot
            main_label=" ",  # Optional main label for the track group
            plddt_label="pLDDT",  # Label for the pLDDT row
            alphamissense_label="AM",  # Label for the AlphaMissense row
            sub_track_height=0.05,  # Height of each row
            sub_track_spacing=0.03,  # Space between rows
        )
        # Plot
        print("\nPlotting...")
        plot_protein_tracks(
            protein_id=uniprot_id,
            sequence_length=seq_length,
            tracks=[
                axis_trk_zoom,
                pdb_trk_detail_zoom,
                ligand_trk,
                alphafold_trk_both,
                ptms_track,
            ],
            view_start_aa=zoom_s,
            view_end_aa=zoom_e,
            figure_width=8,  # May want smaller figure for zoomed view
        )
        print("Done.")

    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
