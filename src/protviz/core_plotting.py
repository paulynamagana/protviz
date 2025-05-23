# src/protviz/core_plotting.py
from typing import List, Optional

import matplotlib.pyplot as plt

from .tracks.base_track import BaseTrack  # Assuming BaseTrack is in this relative path


def plot_protein_tracks(
    protein_id: str,
    sequence_length: int,
    tracks: List[BaseTrack],
    figure_width: float = 10.0,
    figure_height: Optional[float] = None,  # New: Allow direct figure height override
    view_start_aa: Optional[int] = None,
    view_end_aa: Optional[int] = None,
    title_fontsize: int = 12,
    xlabel_fontsize: int = 10,
    save_option: bool = False,
    figure_height_multiplier: float = 3,
):
    """
    Creates a plot with the given protein tracks, supporting zooming.

    Args:
        protein_id: UniProt ID of the protein (for title).
        sequence_length: The total length of the protein sequence.
        tracks: A list of track objects to draw.
        figure_width: Width of the matplotlib figure in inches.
        figure_height: Optional. Exact height of the matplotlib figure in inches.
                       If None, height is calculated based on track content and multiplier.
        view_start_aa: Optional start amino acid for zoomed view. Defaults to 1.
        view_end_aa: Optional end amino acid for zoomed view. Defaults to sequence_length.
        title_fontsize: Font size for the plot title.
        xlabel_fontsize: Font size for the x-axis label.
        save_option: If True, saves the plot to a file instead of showing.
        figure_height_multiplier: Multiplier for total track units to determine figure height,
                                  used only if figure_height is not specified.
    """
    if not tracks:
        print("No tracks provided to plot.")
        return

    # --- View Range Handling ---
    display_start = (
        view_start_aa if view_start_aa is not None and view_start_aa > 0 else 1
    )
    display_end = (
        view_end_aa
        if view_end_aa is not None
        and view_end_aa <= sequence_length
        and view_end_aa > display_start
        else sequence_length
    )
    if sequence_length > 0 and display_start >= display_end:
        display_start = 1
        display_end = sequence_length
    elif sequence_length == 0:
        display_start = 0
        display_end = 1

    if display_start >= display_end and sequence_length > 0:
        print(
            f"Warning: Invalid view range after adjustment. Start AA ({display_start}) must be less than End AA ({display_end}). Plotting full sequence."
        )
        display_start = 1
        display_end = sequence_length
    # --- End View Range Handling ---

    # --- Figure Setup ---
    fig_h: float
    if figure_height is not None:  # User specified exact figure height
        fig_h = figure_height
    else:  # Calculate height based on content and multiplier
        total_plot_height_units = sum(track.get_total_height() for track in tracks)
        calculated_h = total_plot_height_units * figure_height_multiplier

        min_fig_height = 1.0
        if total_plot_height_units <= 0.1:
            fig_h = min_fig_height
        else:
            fig_h = max(min_fig_height, calculated_h)

    fig, ax = plt.subplots(figsize=(figure_width, fig_h))
    # --- End Figure Setup ---

    # --- Track Drawing ---
    current_y_position = 0
    for track in reversed(tracks):
        track.draw(ax, sequence_length, current_y_position, display_start, display_end)
        current_y_position += track.get_total_height()

    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)

    xlim_start = display_start - 0.5
    xlim_end = display_end + 0.5
    if xlim_start >= xlim_end:
        xlim_start = display_start - 1.0
        xlim_end = display_end + 1.0
        if xlim_start < 0 and display_start == 1:
            xlim_start = 0

    ax.set_xlim(xlim_start, xlim_end)

    # --- Labels and Title ---
    ax.set_xlabel("Sequence Position", fontsize=xlabel_fontsize)

    title_str = f"Protein: {protein_id}"
    # Check if a zoom was effectively applied (requested or defaulted to full range)
    is_zoomed_view = (view_start_aa is not None and view_start_aa != 1) or (
        view_end_aa is not None and view_end_aa != sequence_length
    )

    if is_zoomed_view:
        title_str += (
            f" (View: {display_start}-{display_end} aa / Total: {sequence_length} aa)"
        )
    else:
        title_str += f" (Length: {sequence_length} aa)"
    ax.set_title(title_str, fontsize=title_fontsize)
    # --- End Labels and Title ---

    plt.tight_layout(pad=0.8)

    if save_option:
        filename = f"{protein_id.replace('/', '_')}_plot.png"
        plt.savefig(filename, dpi=300)
        print(f"Plot saved as {filename}")
        # plt.close(fig) # Optional: close figure after saving if running in a script

    plt.show()
