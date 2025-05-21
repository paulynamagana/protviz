from protviz import plot_protein_tracks
from protviz.data_retrieval import PDBeClient, TEDClient, get_protein_sequence_length
from protviz.tracks import AxisTrack, LigandInteractionTrack, PDBTrack, TEDDomainsTrack


def main():
    uniprot_id = "O15245"
    pdbe_client = PDBeClient()
    ted_client = TEDClient()

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
        if pdb_coverage:
            print(f"Found {len(pdb_coverage)} PDB coverage entries.")
        else:
            print(f"No PDB coverage data found for {uniprot_id}.")

        # --- Example for TED Domains Track ---
        print(f"\nFetching TED annotations for {uniprot_id}...")
        ted_data = ted_client.get_TED_annotations(uniprot_id)
        if ted_data:
            print(f"Found {len(ted_data)} TED domain annotations.")
            # for ann in ted_data:
            #     print(f"  Chopping: {ann.get('chopping')}, Segments: {ann.get('num_segments')}, CATH: {ann.get('cath_label')}")
        else:
            print(f"No TED annotation data found for {uniprot_id}.")

        # Create tracks
        axis_trk = AxisTrack(
            sequence_length=seq_length,
            label="Sequence",
            tick_interval=max(50, seq_length // 10),
        )

        # Ligand track
        ligand_trk = LigandInteractionTrack(
            interaction_data=ligand_data,
            label="Ligand \n Binding Sites",
            show_ligand_labels=True,  # Show ligand IDs next to their interaction sites
            plotting_option="full",
            site_height=0.2,
            ligand_spacing=0.2,
        )

        # PDB track
        pdb_trk = PDBTrack(
            pdb_data=pdb_coverage,
            label="PDB (Collapsed)",
            plotting_option="collapse",
            bar_height=0.2,
            color="darkturquoise",
        )

        ted_domains_trk = TEDDomainsTrack(
            ted_annotations=ted_data,
            label="TED Domains (Full)",
            plotting_option="full",
            show_domain_labels=True,
        )

        # all tracks
        plot_protein_tracks(
            protein_id=uniprot_id,
            sequence_length=seq_length,
            tracks=[axis_trk, ted_domains_trk, pdb_trk, ligand_trk],
            figure_width=16,
        )

        print("Done with ligand example.")

    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
