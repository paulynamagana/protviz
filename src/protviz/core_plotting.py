# --- src/protviz/core_plotting.py ---
from typing import List, Optional

import matplotlib.pyplot as plt

from .tracks.base_track import BaseTrack


def plot_protein_tracks(
    protein_id: str,
    sequence_length: int,
    tracks: List[BaseTrack],
    figure_width: float = 10.0,
    view_start_aa: Optional[int] = None,
    view_end_aa: Optional[int] = None,
    title_fontsize: int = 12,
    xlabel_fontsize: int = 10,
    save_option: Optional[bool] = False,
):
    """
    Creates a plot with the given protein tracks, supporting zooming.

    Args:
        protein_id: UniProt ID of the protein (for title).
        sequence_length: The total length of the protein sequence.
        tracks: A list of track objects to draw.
        figure_width: Width of the matplotlib figure.
        view_start_aa: Optional start amino acid for zoomed view. Defaults to 1.
        view_end_aa: Optional end amino acid for zoomed view. Defaults to sequence_length.
        title_fontsize: Font size for the plot title.
        xlabel_fontsize: Font size for the x-axis label.
    """
    if not tracks:
        print("No tracks provided to plot.")
        return

    # --- View Range Handling ---
    # Determine the display range for the x-axis
    display_start = (
        view_start_aa if view_start_aa is not None and view_start_aa > 0 else 1
    )
    display_end = (
        view_end_aa
        if view_end_aa is not None and view_end_aa <= sequence_length
        else sequence_length
    )

    if display_start >= display_end:
        print(
            f"Warning: Invalid view range. Start AA ({display_start}) must be less than End AA ({display_end}). Plotting full sequence."
        )
        display_start = 1
        display_end = sequence_length
    # --- End View Range Handling ---

    # --- Figure Setup ---
    # Calculate total height needed for all tracks
    total_plot_height_units = sum(track.get_total_height() for track in tracks)
    figure_height = total_plot_height_units * 5  # Heuristic, adjust as needed
    if figure_height < 2:
        figure_height = 2  # Minimum figure height

    fig, ax = plt.subplots(figsize=(figure_width, figure_height))
    # --- End Figure Setup ---

    # --- Track Drawing ---
    current_y_position = 0  # Tracks are drawn from bottom up
    for track in reversed(tracks):  # Draw from bottom-most track upwards on plot
        # Each track is responsible for drawing its content within the allocated space.
        # Pass the display_start and display_end to each track's draw method
        track.draw(ax, sequence_length, current_y_position, display_start, display_end)
        current_y_position += track.get_total_height()

    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)

    # Set x-axis limits to the determined display range with padding
    ax.set_xlim(
        display_start - 0.5, display_end + 0.5
    )  # Add 0.5 padding for features at exact start/end

    # --- Labels and Title ---
    ax.set_xlabel("Sequence Position", fontsize=xlabel_fontsize)

    title_str = f"Protein: {protein_id}"
    if view_start_aa is not None or view_end_aa is not None:
        # Indicate zoomed view in the title
        title_str += (
            f" (View: {display_start}-{display_end} aa / Total: {sequence_length} aa)"
        )
    else:
        # Indicate full sequence view in the title
        title_str += f" (Length: {sequence_length} aa)"
    ax.set_title(title_str, fontsize=title_fontsize)
    # --- End Labels and Title ---

    plt.tight_layout(pad=0.5)  # Adjust layout to prevent labels from overlapping
    if save_option is True:
        plt.savefig(f"{protein_id}_plot.png", dpi=300)
        print(f"Plot saved as {protein_id}_plot.png")
        plt.show()
    else:
        plt.show()
