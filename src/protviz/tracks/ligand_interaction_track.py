# --- src/protviz/tracks/ligand_interaction_track.py
from typing import Any, Dict, List, Literal

import matplotlib.patches as patches
import matplotlib.pyplot as plt

from ..utils.styling import (
    ANNOTATION_TEXT_SIZE,
    TRACK_LABEL_COLOR,
    TRACK_LABEL_FONT_SIZE,
    get_distinct_colors,
)
from .base_track import BaseTrack

PlottingOptionLigand = Literal["full", "collapse"]


class LigandInteractionTrack(BaseTrack):
    def __init__(
        self,
        interaction_data: List[Dict[str, Any]],
        label: str = "Ligand Binding",
        plotting_option: PlottingOptionLigand = "full",
        site_height: float = 0.1,
        ligand_spacing: float = 0.1,
        default_color: str = "purple",
        show_ligand_labels: bool = True,
        ligand_label_offset: float = 0.015,
        full_mode_padding: float = 0.1,  # Padding for "full" mode
        collapse_mode_padding: float = 0.1,
    ):  # Smaller padding for "collapse"
        super().__init__(label=label, height=0, padding=0.1)

        self.interaction_data_orig = interaction_data
        self.plotting_option = plotting_option
        self.site_height = site_height  # This is content height for collapse, or sub-bar height for full
        self.ligand_spacing = ligand_spacing
        self.default_color = default_color
        self.show_ligand_labels = show_ligand_labels
        self.ligand_label_offset = ligand_label_offset

        # Store mode-specific paddings
        self.full_mode_padding = full_mode_padding
        self.collapse_mode_padding = collapse_mode_padding

        # Set initial padding based on mode - this will be used by get_total_height() via BaseTrack
        if self.plotting_option == "collapse":
            self.padding = self.collapse_mode_padding
        else:  # full mode
            self.padding = self.full_mode_padding

        self._grouped_by_ligand_view: Dict[str, List[Dict[str, Any]]] = {}
        self._unique_ligands_view: List[str] = []
        self._num_ligand_lanes_view: int = 0
        self._ligand_colors_view: Dict[str, str] = {}
        self._merged_ligand_sites_view: List[Dict[str, int]] = []

        if not self.interaction_data_orig:
            # If no data, content height is effectively site_height for message display
            self.height = self.site_height

    # _extract_all_sites and _calculate_merged_site_regions remain the same
    def _extract_all_sites(
        self, interaction_data: List[Dict[str, Any]]
    ) -> List[Dict[str, int]]:
        all_sites = []
        for interaction_context in interaction_data:
            for site_segment in interaction_context.get(
                "binding_site_uniprot_residues", []
            ):
                if (
                    site_segment.get("startIndex") is not None
                    and site_segment.get("endIndex") is not None
                ):
                    all_sites.append(
                        {
                            "startIndex": site_segment["startIndex"],
                            "endIndex": site_segment["endIndex"],
                        }
                    )
        return sorted(all_sites, key=lambda x: (x["startIndex"], x["endIndex"]))

    def _calculate_merged_site_regions(
        self, all_sites: List[Dict[str, int]]
    ) -> List[Dict[str, int]]:
        if not all_sites:
            return []
        merged = []
        current_start, current_end = -1, -1
        for site in all_sites:
            site_start, site_end = site["startIndex"], site["endIndex"]
            if current_start == -1:
                current_start, current_end = site_start, site_end
            else:
                if site_start <= current_end + 1:
                    current_end = max(current_end, site_end)
                else:
                    merged.append(
                        {"startIndex": current_start, "endIndex": current_end}
                    )
                    current_start, current_end = site_start, site_end
        if current_start != -1:
            merged.append({"startIndex": current_start, "endIndex": current_end})
        return merged

    def _layout_entries_for_view(self, view_start_aa: int, view_end_aa: int):
        """Calculates layout and sets self.height and self.padding based on plotting_option and current view."""

        # Set padding based on mode first
        if self.plotting_option == "collapse":
            self.padding = self.collapse_mode_padding
        else:  # full mode
            self.padding = self.full_mode_padding

        temp_grouped_by_ligand: Dict[str, List[Dict[str, Any]]] = {}
        for interaction_context in self.interaction_data_orig:
            ligand_id = interaction_context.get("ligand_id", "UNKNOWN_LIGAND")
            if ligand_id not in temp_grouped_by_ligand:
                temp_grouped_by_ligand[ligand_id] = []
            for site_segment in interaction_context.get(
                "binding_site_uniprot_residues", []
            ):
                s, e = site_segment.get("startIndex"), site_segment.get("endIndex")
                if (
                    s is not None
                    and e is not None
                    and not (e < view_start_aa or s > view_end_aa)
                ):
                    temp_grouped_by_ligand[ligand_id].append(site_segment)

        self._grouped_by_ligand_view = {
            lig_id: sites for lig_id, sites in temp_grouped_by_ligand.items() if sites
        }
        self._unique_ligands_view = sorted(list(self._grouped_by_ligand_view.keys()))
        self._num_ligand_lanes_view = len(self._unique_ligands_view)
        self._ligand_colors_view = {
            lig_id: color
            for lig_id, color in zip(
                self._unique_ligands_view,
                get_distinct_colors(
                    self._num_ligand_lanes_view, self._num_ligand_lanes_view
                ),
            )
        }

        if self.plotting_option == "full":
            if self._num_ligand_lanes_view > 0:
                # Content height for full mode
                self.height = (self._num_ligand_lanes_view * self.site_height) + (
                    max(0, self._num_ligand_lanes_view - 1) * self.ligand_spacing
                )
            else:
                self.height = self.site_height  # For "No data" message

        elif self.plotting_option == "collapse":
            all_visible_sites = []
            for ligand_id in self._unique_ligands_view:
                for site_segment in self._grouped_by_ligand_view[ligand_id]:
                    all_visible_sites.append(site_segment)

            all_visible_sites_sorted = sorted(
                all_visible_sites, key=lambda x: (x["startIndex"], x["endIndex"])
            )
            self._merged_ligand_sites_view = self._calculate_merged_site_regions(
                all_visible_sites_sorted
            )
            # Content height for collapse mode is just the site_height
            self.height = self.site_height
        else:
            raise ValueError(f"Unknown plotting_option: {self.plotting_option}")

        # Ensure minimum height if no actual data to draw (for the "No data in view" message)
        if (
            not self.interaction_data_orig
            or (self.plotting_option == "full" and self._num_ligand_lanes_view == 0)
            or (
                self.plotting_option == "collapse"
                and not self._merged_ligand_sites_view
            )
        ):
            self.height = self.site_height

    def draw(
        self,
        ax: plt.Axes,
        protein_sequence_length: int,
        track_y_bottom: float,
        view_start_aa: int,
        view_end_aa: int,
    ):
        self._layout_entries_for_view(view_start_aa, view_end_aa)

        # Note: self.padding is now set by _layout_entries_for_view
        content_block_y_start = track_y_bottom + self.padding

        no_data_to_draw = False
        if self.plotting_option == "full" and not self._unique_ligands_view:
            no_data_to_draw = True
        elif self.plotting_option == "collapse" and not self._merged_ligand_sites_view:
            no_data_to_draw = True
        if not self.interaction_data_orig:
            no_data_to_draw = True

        if no_data_to_draw:
            message_y_center = content_block_y_start + self.height / 2
            ax.text(
                (view_start_aa + view_end_aa) / 2,
                message_y_center,
                "No ligand data in view",
                ha="center",
                va="center",
                fontsize=ANNOTATION_TEXT_SIZE - 1,
                color="gray",
                style="italic",
            )
            if self.label:
                ax.text(
                    0,
                    message_y_center,
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

        track_content_center_y = content_block_y_start + self.height / 2
        if self.label:
            ax.text(
                0,
                track_content_center_y,
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
        for i, ligand_id in enumerate(self._unique_ligands_view):
            lane_index = i
            lane_y_start = content_block_y_start + lane_index * (
                self.site_height + self.ligand_spacing
            )

            ligand_color = self._ligand_colors_view.get(ligand_id, self.default_color)
            interaction_sites = self._grouped_by_ligand_view.get(ligand_id, [])

            has_visible_sites_for_this_ligand = False
            for site in interaction_sites:
                site_start = site.get("startIndex")
                site_end = site.get("endIndex")

                plot_s = max(site_start, view_start_aa)
                plot_e = min(site_end, view_end_aa)

                if plot_s > plot_e:
                    continue
                has_visible_sites_for_this_ligand = True

                bar_width = plot_e - plot_s + 1
                if bar_width <= 0:
                    bar_width = 0.8
                rect_x = (
                    plot_s - 0.4
                    if bar_width == 0.8 and (plot_e - plot_s + 1) == 1
                    else plot_s
                )

                rect = patches.Rectangle(
                    (rect_x, lane_y_start),
                    bar_width,
                    self.site_height,
                    linewidth=0.5,
                    edgecolor="black",
                    facecolor=ligand_color,
                    alpha=0.8,
                    clip_on=True,
                )
                ax.add_patch(rect)

            if self.show_ligand_labels and has_visible_sites_for_this_ligand:
                label_y_center = lane_y_start + self.site_height / 2
                ax.text(
                    1.0 + self.ligand_label_offset,
                    label_y_center,
                    ligand_id,
                    ha="left",
                    va="center",
                    fontsize=ANNOTATION_TEXT_SIZE,
                    color=ligand_color,
                    transform=ax.get_yaxis_transform(),
                    clip_on=False,
                )

    def _draw_collapse_mode(
        self,
        ax: plt.Axes,
        view_start_aa: int,
        view_end_aa: int,
        content_block_y_start: float,
    ):
        lane_y_start = content_block_y_start

        for region in self._merged_ligand_sites_view:
            region_start = region["startIndex"]
            region_end = region["endIndex"]

            plot_s = max(region_start, view_start_aa)
            plot_e = min(region_end, view_end_aa)

            if plot_s > plot_e:
                continue

            bar_width = plot_e - plot_s + 1
            if bar_width <= 0:
                bar_width = 0.8
            rect_x = (
                plot_s - 0.4
                if bar_width == 0.8 and (plot_e - plot_s + 1) == 1
                else plot_s
            )

            rect = patches.Rectangle(
                (rect_x, lane_y_start),
                bar_width,
                self.site_height,
                linewidth=0.5,
                edgecolor="black",
                facecolor=self.default_color,
                alpha=0.85,
                clip_on=True,
            )
            ax.add_patch(rect)
