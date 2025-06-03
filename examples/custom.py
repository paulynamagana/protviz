from protviz import plot_protein_tracks
from protviz.data_retrieval import get_protein_sequence_length
from protviz.tracks import AxisTrack, CustomTrack


def main():
    uniprot_id = "Q92794"  # Replace with your protein ID

    seq_length = get_protein_sequence_length(uniprot_id)
    # 1. Define your custom annotation data
    my_feature_data = [
        {
            "position": 25,
            "label": "Active Site Residue",
            "row_label": "Functional",
            "color": "darkgreen",
            "display_type": "marker",
            "marker_symbol": "*",
            "marker_size": 10,
        },
        {
            "start": 50,
            "end": 75,
            "label": "Motif A",
            "row_label": "Sequence Motifs",
            "color": "navy",
        },
        {
            "position": 110,
            "label": "PTM X",
            "row_label": "Modifications",
            "color": "saddlebrown",
            "display_type": "marker",
            "marker_symbol": "P",
        },
        {
            "start": 130,
            "end": 160,
            "label": "Region B",
            "row_label": "Structural Features",
            "color": "teal",
            "display_type": "bar",
        },
        {
            "position": 180,
            "label": "PTM Y",
            "row_label": "Modifications",
            "color": "darkorange",
            "display_type": "marker",
            "marker_symbol": "X",
        },
    ]

    # 2. Create the AxisTrack
    axis_trk = AxisTrack(sequence_length=seq_length, label="Sequence")

    # 3. Create the CustomTrack
    custom_display_track = CustomTrack(
        annotation_data=my_feature_data,
        show_row_labels=True,
        show_ann_labels=True,
        ann_height=0.05,
        ann_spacing=0.08,
    )

    # 4. Plot
    plot_protein_tracks(
        protein_id=uniprot_id,
        sequence_length=seq_length,
        tracks=[axis_trk, custom_display_track],
        figure_width=10,
        figure_height=8,
        save_option=True,
    )


if __name__ == "__main__":
    main()
