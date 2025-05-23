from protviz import plot_protein_tracks
from protviz.data_retrieval import (
    InterProClient,
    PDBeClient,
    get_protein_sequence_length,
)
from protviz.tracks import (
    AxisTrack,
    CustomTrack,
    InterProTrack,
    LigandInteractionTrack,
    PDBTrack,
)


def main():
    uniprot_id = "Q92794"
    pdbe_client = PDBeClient()
    interpro_client = InterProClient()

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

        print(f"Fetching Pfam annotations for {uniprot_id}...")
        pfam_domain_data = interpro_client.get_pfam_annotations(uniprot_id)

        print(f"Fetching CATH-Gene3D annotations for {uniprot_id}...")
        cath_domain_data = interpro_client.get_cathgene3d_annotations(uniprot_id)

        combined_data = [
            {
                "start": 12,
                "end": 15,
                "label": "KRAS G12V",
                "row_label": "Mutations",
                "color": "red",
            },
            {
                "start": 50,
                "end": 150,
                "label": "Ras Domain",
                "row_label": "New domain",
                "color": "blue",
            },
            {
                "position": 18,
                "label": "Drug Binding",
                "row_label": "Interactions",
                "color": "green",
            },
        ]

        pfam_track_full = InterProTrack(
            domain_data=pfam_domain_data,
            database_name_for_label="Pfam",
            label="Pfam (Full)",
            plotting_option="full",
        )

        cath_track_full = InterProTrack(
            domain_data=cath_domain_data,
            database_name_for_label="CATH-Gene3D",
            label="CATH-Gene3D (Full)",
            plotting_option="full",
        )

        ptms_track = CustomTrack(
            annotation_data=combined_data,
            label="",  # Main track label
            ann_height=0.3,  # Height for the single row of PTMs
            show_row_labels=True,
            padding=0.05,
        )

        # Create tracks
        axis_trk = AxisTrack(
            sequence_length=seq_length,
            label="Sequence",
        )

        ligand_trk = LigandInteractionTrack(
            interaction_data=ligand_data,
            label="Ligand Binding Sites",
            show_ligand_labels=True,  # Show ligand IDs next to their interaction sites
            plotting_option="collapse",
            site_height=0.2,
            ligand_spacing=0.2,
        )

        pdb_trk_lig = PDBTrack(
            pdb_data=pdb_coverage,
            label="PDB (Collapsed)",
            plotting_option="collapse",
            bar_height=0.2,
            color="darkturquoise",
        )

        plot_protein_tracks(
            protein_id=uniprot_id,
            sequence_length=seq_length,
            tracks=[
                axis_trk,
                pdb_trk_lig,
                ligand_trk,
                cath_track_full,
                pfam_track_full,
                ptms_track,
            ],
            figure_width=12,
            figure_height=5,
            save_option=True,
        )

    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
