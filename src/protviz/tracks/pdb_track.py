# --- src/protviz/tracks/pdb_track.py ---
from typing import Any, Dict, List, Literal

import matplotlib.patches as patches
import matplotlib.pyplot as plt

from ..utils.styling import (
    ANNOTATION_TEXT_SIZE,
    TRACK_LABEL_COLOR,
    TRACK_LABEL_FONT_SIZE,
)
from .base_track import BaseTrack

PlottingOption = Literal["full", "collapse"]

PlottingOptionPDB = Literal["full", "collapse"]


class PDBTrack(BaseTrack):
    def __init__(
        self,
        pdb_data: List[Dict[str, Any]],
        label: str = "PDB Coverage",
        plotting_option: PlottingOptionPDB = "full",
        bar_height: float = 0.1,
        bar_vertical_spacing: float = 0.05,
        color: str = "skyblue",
        text_color: str = "black",
        text_size: int = ANNOTATION_TEXT_SIZE,  # Use from styling
        full_mode_padding: float = 0.1,  # Padding for "full" mode
        collapse_mode_padding: float = 0.1,  # Smaller padding for "collapse"
    ):
        # Initialize with temporary padding, will be set correctly in _layout_entries_for_view
        super().__init__(
            label=label,
            height=0,
            padding=full_mode_padding
            if plotting_option == "full"
            else collapse_mode_padding,
        )

        self.pdb_data_orig = sorted(
            pdb_data, key=lambda x: (x.get("unp_start", 0), x.get("unp_end", 0))
        )
        self.plotting_option = plotting_option
        self.bar_height = (
            bar_height  # Content height for collapse, sub-bar height for full
        )
        self.bar_vertical_spacing = bar_vertical_spacing
        self.color = color
        self.text_color = text_color
        self.text_size = text_size

        self.full_mode_padding = full_mode_padding
        self.collapse_mode_padding = collapse_mode_padding

        self._lanes_assignment_view: List[tuple[Dict[str, Any], int]] = []
        self._merged_regions_view: List[Dict[str, int]] = []
        self._num_lanes_view: int = 0

        if not self.pdb_data_orig:
            self.height = self.bar_height

    def _calculate_merged_regions(
        self, current_pdb_data: List[Dict[str, Any]]
    ) -> List[Dict[str, int]]:
        if not current_pdb_data:
            return []
        pdb_data_sorted = sorted(
            current_pdb_data, key=lambda x: (x.get("unp_start", 0), x.get("unp_end", 0))
        )
        merged = []
        if not pdb_data_sorted:
            return merged

        # Handle case where pdb_data_sorted might be empty after filtering, though caught by outer if
        if not pdb_data_sorted:
            return []

        current_start, current_end = (
            pdb_data_sorted[0].get("unp_start", 0),
            pdb_data_sorted[0].get("unp_end", 0),
        )
        # Ensure current_start and current_end are valid before proceeding
        if (
            current_start == 0 and current_end == 0 and len(pdb_data_sorted) == 1
        ):  # Only one invalid entry
            return []
        if (
            current_start == 0 and current_end == 0
        ):  # If first entry is invalid, try to find first valid
            found_valid_start = False
            for entry_init in pdb_data_sorted:
                cs_init, ce_init = (
                    entry_init.get("unp_start", 0),
                    entry_init.get("unp_end", 0),
                )
                if not (cs_init == 0 and ce_init == 0):
                    current_start, current_end = cs_init, ce_init
                    found_valid_start = True
                    break
            if not found_valid_start:
                return []  # No valid entries at all

        for i in range(1, len(pdb_data_sorted)):
            next_start, next_end = (
                pdb_data_sorted[i].get("unp_start", 0),
                pdb_data_sorted[i].get("unp_end", 0),
            )
            if next_start == 0 and next_end == 0:
                continue  # Skip invalid entries

            if next_start <= current_end + 1:
                current_end = max(current_end, next_end)
            else:
                merged.append({"unp_start": current_start, "unp_end": current_end})
                current_start, current_end = next_start, next_end
        merged.append({"unp_start": current_start, "unp_end": current_end})
        return merged

    def _layout_entries_for_view(self, view_start_aa: int, view_end_aa: int):
        if self.plotting_option == "collapse":
            self.padding = self.collapse_mode_padding
        else:
            self.padding = self.full_mode_padding

        visible_pdb_data = [
            e
            for e in self.pdb_data_orig
            if not (
                e.get("unp_end", 0) < view_start_aa
                or e.get("unp_start", 0) > view_end_aa
            )
        ]

        if not visible_pdb_data:
            self._num_lanes_view = 0  # No lanes if no visible data
            self.height = self.bar_height  # For "No data" message
            self._lanes_assignment_view = []
            self._merged_regions_view = []
            return

        if self.plotting_option == "full":
            lanes_end_pos: List[int] = []
            self._lanes_assignment_view = []
            # Sort visible_pdb_data for consistent lane assignment if order matters
            # visible_pdb_data_sorted = sorted(visible_pdb_data, key=lambda x: (x.get('unp_start', 0), x.get('unp_end', 0)))
            for entry in (
                visible_pdb_data
            ):  # Use visible_pdb_data directly as it's already filtered
                layout_s, layout_e = (
                    max(entry.get("unp_start", 0), view_start_aa),
                    min(entry.get("unp_end", 0), view_end_aa),
                )
                if layout_s >= layout_e:
                    continue
                placed = False
                for i, lane_end_coord in enumerate(lanes_end_pos):
                    if layout_s > lane_end_coord:
                        self._lanes_assignment_view.append((entry, i))
                        lanes_end_pos[i] = layout_e
                        placed = True
                        break
                if not placed:
                    self._lanes_assignment_view.append((entry, len(lanes_end_pos)))
                    lanes_end_pos.append(layout_e)

            self._num_lanes_view = (
                len(lanes_end_pos) if lanes_end_pos else 0
            )  # No lanes if no assignments
            if self._num_lanes_view > 0:
                self.height = (self._num_lanes_view * self.bar_height) + (
                    max(0, self._num_lanes_view - 1) * self.bar_vertical_spacing
                )
            else:  # No items were placed in lanes (e.g., all had zero visible width)
                self.height = self.bar_height  # For "No data" message

        elif self.plotting_option == "collapse":
            self._merged_regions_view = self._calculate_merged_regions(visible_pdb_data)
            self._num_lanes_view = 1  # Conceptually one lane for collapse mode
            self.height = self.bar_height  # Content height is just the bar height

        # If after layout, no actual elements to draw, still ensure height for message
        if (self.plotting_option == "full" and not self._lanes_assignment_view) or (
            self.plotting_option == "collapse" and not self._merged_regions_view
        ):
            self.height = self.bar_height

    def draw(
        self,
        ax: plt.Axes,
        protein_sequence_length: int,
        track_y_bottom: float,
        view_start_aa: int,
        view_end_aa: int,
    ):
        self._layout_entries_for_view(view_start_aa, view_end_aa)
        content_block_y_start = track_y_bottom + self.padding

        no_data = False
        if self.plotting_option == "full" and not self._lanes_assignment_view:
            no_data = True
        if self.plotting_option == "collapse" and not self._merged_regions_view:
            no_data = True
        if not self.pdb_data_orig:
            no_data = True  # If original data was empty

        if no_data:
            msg_y = content_block_y_start + self.height / 2
            ax.text(
                (view_start_aa + view_end_aa) / 2,
                msg_y,
                "No PDB data in view",
                ha="center",
                va="center",
                fontsize=ANNOTATION_TEXT_SIZE - 1,
                color="gray",
                style="italic",
            )
            if self.label:
                ax.text(
                    0,
                    msg_y,
                    self.label,
                    ha="right",
                    va="center",
                    fontsize=TRACK_LABEL_FONT_SIZE,
                    color=TRACK_LABEL_COLOR,
                    transform=ax.get_yaxis_transform(),
                    clip_on=False,
                )
            return

        if self.plotting_option == "full":
            self._draw_full_mode(ax, view_start_aa, view_end_aa, content_block_y_start)
        elif self.plotting_option == "collapse":
            self._draw_collapse_mode(
                ax, view_start_aa, view_end_aa, content_block_y_start
            )

        if self.label:
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

    def _draw_full_mode(
        self,
        ax: plt.Axes,
        view_start_aa: int,
        view_end_aa: int,
        content_block_y_start: float,
    ):
        for entry, lane_idx in self._lanes_assignment_view:
            plot_s, plot_e = (
                max(entry.get("unp_start", 0), view_start_aa),
                min(entry.get("unp_end", 0), view_end_aa),
            )
            if plot_s >= plot_e:
                continue
            bar_w = plot_e - plot_s
            if bar_w <= 0 and plot_e == plot_s:
                bar_w = 0.8  # Min width for single residue

            bar_y = content_block_y_start + lane_idx * (
                self.bar_height + self.bar_vertical_spacing
            )
            rect = patches.Rectangle(
                (plot_s, bar_y),
                bar_w,
                self.bar_height,
                linewidth=0.5,
                edgecolor=self.text_color,
                facecolor=self.color,
                alpha=0.75,
                clip_on=True,
            )
            ax.add_patch(rect)
            pdb_id = entry.get("pdb_id", "N/A")
            text_y = bar_y + self.bar_height / 2
            view_span = view_end_aa - view_start_aa + 1
            min_bar_w_text = len(pdb_id) * 0.015 * view_span  # Heuristic
            if bar_w > min_bar_w_text:
                ax.text(
                    plot_s + bar_w / 2,
                    text_y,
                    pdb_id,
                    ha="center",
                    va="center",
                    fontsize=self.text_size,
                    color=self.text_color,
                    clip_on=True,
                )
            elif (
                bar_w > 0.001 * view_span and view_span > 0
            ):  # Avoid division by zero if view_span is 0
                label_offset_x = (
                    view_span * 0.005 if view_span > 0 else 0.5
                )  # Default offset if view_span is 0
                text_x = plot_e + label_offset_x
                ha_ext = "left"
                # Estimate text width to avoid overlap (very rough)
                estimated_text_width = (
                    len(pdb_id)
                    * (self.text_size / 72)
                    * (ax.get_figure().get_dpi())
                    * (
                        view_span / ax.get_xlim()[1] - ax.get_xlim()[0]
                        if ax.get_xlim()[1] > ax.get_xlim()[0]
                        else 1
                    )
                )

                if (
                    text_x + estimated_text_width > view_end_aa
                    and plot_s - label_offset_x - estimated_text_width > view_start_aa
                ):  # check if it fits on left
                    text_x, ha_ext = plot_s - label_offset_x, "right"
                elif (
                    text_x + estimated_text_width > view_end_aa
                ):  # If still doesn't fit left, don't draw or try other strategy
                    continue  # Or try to place above/below, or shorten text

                if not (text_x < view_start_aa and ha_ext == "right"):
                    ax.text(
                        text_x,
                        text_y,
                        pdb_id,
                        ha=ha_ext,
                        va="center",
                        fontsize=max(1, self.text_size - 1),
                        color=self.text_color,
                        clip_on=True,
                    )

    def _draw_collapse_mode(
        self,
        ax: plt.Axes,
        view_start_aa: int,
        view_end_aa: int,
        content_block_y_start: float,
    ):
        bar_y = content_block_y_start
        for region in self._merged_regions_view:
            plot_s, plot_e = (
                max(region.get("unp_start", 0), view_start_aa),
                min(region.get("unp_end", 0), view_end_aa),
            )
            if plot_s >= plot_e:
                continue
            bar_w = plot_e - plot_s
            if bar_w <= 0 and plot_e == plot_s:
                bar_w = 0.8
            rect = patches.Rectangle(
                (plot_s, bar_y),
                bar_w,
                self.bar_height,
                linewidth=0.5,
                edgecolor=self.text_color,
                facecolor=self.color,
                alpha=0.85,
                clip_on=True,
            )
            ax.add_patch(rect)
