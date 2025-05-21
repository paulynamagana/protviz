# --- src/protviz/tracks/alphafold_track.py ---
import logging
from collections import defaultdict
from typing import Any, Dict, List, Optional

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from ..utils.styling import (
    ANNOTATION_TEXT_SIZE,
    TRACK_LABEL_COLOR,
    TRACK_LABEL_FONT_SIZE,
    get_alphamissense_color,
    get_plddt_color,
)
from .base_track import BaseTrack

logger_af_track = logging.getLogger(__name__)


class AlphaFoldTrack(BaseTrack):
    """
    A track to display metrics from AlphaFold predictions, such as pLDDT
    and average AlphaMissense pathogenicity scores.
    """

    def __init__(
        self,
        afdb_data: Dict[
            str, List[Dict[str, Any]]
        ],  # Expects {"plddt": [...], "alphamissense": [...]}
        plotting_options: Optional[
            List[str]
        ] = None,  # e.g., ["plddt", "alphamissense"]
        plddt_label: str = "pLDDT",
        alphamissense_label: str = "AlphaMissense (avg)",
        sub_track_height: float = 0.4,  # Height for each individual data type row
        sub_track_spacing: float = 0.1,  # Vertical spacing between sub-tracks
        main_label: Optional[str] = "AlphaFold Data",
    ):  # Overall label for the track group
        """
        Initializes the AlphaFoldTrack.

        Args:
            afdb_data: Dictionary containing AlphaFold data lists, keyed by data type
                       (e.g., "plddt", "alphamissense").
            plotting_options: A list of strings specifying which data types to plot.
                              If None, attempts to plot all available data in afdb_data.
            plddt_label: Label for the pLDDT sub-track.
            alphamissense_label: Label for the AlphaMissense sub-track.
            sub_track_height: Height for each row (pLDDT, AlphaMissense).
            sub_track_spacing: Vertical space between rows if multiple are plotted.
            main_label: Optional overall label for the group of AlphaFold tracks.
        """

        self.afdb_data = afdb_data if afdb_data else {}
        self.plddt_data = self.afdb_data.get("plddt", [])
        self.alphamissense_data_raw = self.afdb_data.get("alphamissense", [])

        if plotting_options is None:
            self.plotting_options = [
                opt for opt in ["plddt", "alphamissense"] if self.afdb_data.get(opt)
            ]
        else:
            self.plotting_options = [
                opt for opt in plotting_options if self.afdb_data.get(opt)
            ]

        self.plddt_label = plddt_label
        self.alphamissense_label = alphamissense_label
        self.sub_track_height = sub_track_height
        self.sub_track_spacing = sub_track_spacing

        # Calculate total height based on active sub-tracks
        num_active_sub_tracks = len(self.plotting_options)
        calculated_content_height = 0
        if num_active_sub_tracks > 0:
            calculated_content_height = (
                num_active_sub_tracks * self.sub_track_height
            ) + (max(0, num_active_sub_tracks - 1) * self.sub_track_spacing)
        else:  # No valid plotting options or no data
            calculated_content_height = self.sub_track_height  # For "No data" message

        super().__init__(
            label=main_label, height=calculated_content_height
        )  # main_label is for the whole block

        # Pre-calculate average AlphaMissense scores if needed
        self.avg_am_scores_per_residue: Dict[int, float] = {}
        if "alphamissense" in self.plotting_options and self.alphamissense_data_raw:
            self.avg_am_scores_per_residue = (
                self._calculate_average_am_pathogenicity_per_residue(
                    self.alphamissense_data_raw
                )
            )

        if not self.plotting_options:
            logger_af_track.warning(
                "No valid AlphaFold data types to plot or provided data is empty."
            )

    def _calculate_average_am_pathogenicity_per_residue(
        self, am_data_list: List[Dict[str, Any]]
    ) -> Dict[int, float]:
        """
        Calculates the average AlphaMissense pathogenicity score for each residue.
        """
        if not am_data_list:
            logger_af_track.info("No AlphaMissense data provided for averaging.")
            return {}

        scores_by_residue = defaultdict(list)
        for entry in am_data_list:
            try:
                residue = entry["residue_number"]
                score = entry["am_pathogenicity"]
                scores_by_residue[residue].append(score)
            except KeyError as e:
                logger_af_track.warning(
                    f"Missing key in AlphaMissense entry: {e} - Entry: {entry}"
                )
                continue
            except TypeError as e:  # Handle if score is not a number
                logger_af_track.warning(
                    f"Invalid score type for AlphaMissense entry: {e} - Entry: {entry}"
                )
                continue

        avg_scores = {}
        for residue, scores in scores_by_residue.items():
            if scores:  # Ensure scores list is not empty
                avg_scores[residue] = np.mean(scores)

        if not avg_scores:
            logger_af_track.info(
                "No valid scores found to calculate average AlphaMissense pathogenicity."
            )
        return avg_scores

    def _draw_plddt_sub_track(
        self,
        ax: plt.Axes,
        y_start_sub_track: float,
        view_start_aa: int,
        view_end_aa: int,
    ):
        """Draws the pLDDT scores as a colored bar."""
        if not self.plddt_data:
            ax.text(
                (view_start_aa + view_end_aa) / 2,
                y_start_sub_track + self.sub_track_height / 2,
                "No pLDDT data",
                ha="center",
                va="center",
                fontsize=ANNOTATION_TEXT_SIZE - 1,
                color="gray",
            )
            return

        # Create a bar for each residue in the view
        for score_entry in self.plddt_data:
            res_num = score_entry.get("residue_number")
            plddt = score_entry.get("plddt")

            if res_num is None or plddt is None:
                continue

            # Plot only if residue is within the current view
            if view_start_aa <= res_num <= view_end_aa:
                color = get_plddt_color(plddt)
                # Draw a small rectangle for each residue. x is res_num, width is 1.
                # Adjust x to be start of residue, not center, if width is 1.
                rect_x = res_num - 0.5  # Center the bar on the residue number
                rect = patches.Rectangle(
                    (rect_x, y_start_sub_track),
                    1,
                    self.sub_track_height,
                    facecolor=color,
                    edgecolor=color,  # No visible edge or slight darker
                    clip_on=True,
                )
                ax.add_patch(rect)

        # Add sub-track label
        ax.text(
            0,
            y_start_sub_track + self.sub_track_height / 2,
            self.plddt_label,
            ha="right",
            va="center",
            fontsize=TRACK_LABEL_FONT_SIZE - 1,
            color=TRACK_LABEL_COLOR,
            transform=ax.get_yaxis_transform(),
            clip_on=False,
        )

    def _draw_alphamissense_sub_track(
        self,
        ax: plt.Axes,
        y_start_sub_track: float,
        view_start_aa: int,
        view_end_aa: int,
    ):
        """Draws the average AlphaMissense pathogenicity scores."""
        if not self.avg_am_scores_per_residue:
            ax.text(
                (view_start_aa + view_end_aa) / 2,
                y_start_sub_track + self.sub_track_height / 2,
                "No AlphaMissense data",
                ha="center",
                va="center",
                fontsize=ANNOTATION_TEXT_SIZE - 1,
                color="gray",
            )
            return

        for res_num in range(view_start_aa, view_end_aa + 1):
            avg_score = self.avg_am_scores_per_residue.get(res_num)
            if avg_score is not None:
                color = get_alphamissense_color(avg_score)
                rect_x = res_num - 0.5
                rect = patches.Rectangle(
                    (rect_x, y_start_sub_track),
                    1,
                    self.sub_track_height,
                    facecolor=color,
                    edgecolor=color,
                    clip_on=True,
                )
                ax.add_patch(rect)

        # Add sub-track label
        ax.text(
            0,
            y_start_sub_track + self.sub_track_height / 2,
            self.alphamissense_label,
            ha="right",
            va="center",
            fontsize=TRACK_LABEL_FONT_SIZE - 1,
            color=TRACK_LABEL_COLOR,
            transform=ax.get_yaxis_transform(),
            clip_on=False,
        )

    def draw(
        self,
        ax: plt.Axes,
        protein_sequence_length: int,
        track_y_bottom: float,  # y-coordinate for the bottom of the entire AlphaFoldTrack block
        view_start_aa: int,
        view_end_aa: int,
    ):
        """Draws the selected AlphaFold data sub-tracks."""

        content_block_y_start = (
            track_y_bottom + self.padding
        )  # Bottom of the drawable area for this track group

        if not self.plotting_options:
            # Display a general "No AlphaFold data" message if no valid options
            ax.text(
                (view_start_aa + view_end_aa) / 2,
                content_block_y_start + self.height / 2,  # Use calculated total height
                "No AlphaFold data to display",
                ha="center",
                va="center",
                fontsize=ANNOTATION_TEXT_SIZE - 1,
                color="gray",
                style="italic",
            )
            if self.label:  # Main track label
                ax.text(
                    0,
                    content_block_y_start + self.height / 2,
                    self.label,
                    ha="right",
                    va="center",
                    fontsize=TRACK_LABEL_FONT_SIZE,
                    color=TRACK_LABEL_COLOR,
                    transform=ax.get_yaxis_transform(),
                    clip_on=False,
                )
            return

        num_drawn_sub_tracks = 0
        y_offset_for_current_sub_track = 0

        if "plddt" in self.plotting_options:
            actual_sub_track_y_start = (
                content_block_y_start + y_offset_for_current_sub_track
            )
            self._draw_plddt_sub_track(
                ax, actual_sub_track_y_start, view_start_aa, view_end_aa
            )
            y_offset_for_current_sub_track += (
                self.sub_track_height + self.sub_track_spacing
            )
            num_drawn_sub_tracks += 1

        if "alphamissense" in self.plotting_options:
            actual_sub_track_y_start = (
                content_block_y_start + y_offset_for_current_sub_track
            )
            self._draw_alphamissense_sub_track(
                ax, actual_sub_track_y_start, view_start_aa, view_end_aa
            )
            # y_offset_for_current_sub_track += self.sub_track_height + self.sub_track_spacing # Not needed after last
            num_drawn_sub_tracks += 1

        # Add the main track label if provided and if any sub-tracks were drawn
        if self.label and num_drawn_sub_tracks > 0:
            main_label_y_center = (
                content_block_y_start + self.height / 2
            )  # Center of the whole block
            ax.text(
                0,
                main_label_y_center,
                self.label,
                ha="right",
                va="center",
                fontsize=TRACK_LABEL_FONT_SIZE,
                color=TRACK_LABEL_COLOR,
                transform=ax.get_yaxis_transform(),
                clip_on=False,
            )
