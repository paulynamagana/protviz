# --- src/protviz/tracks/ted_track.py ---
import logging
from collections import defaultdict
from typing import Any, Dict, List, Literal, Tuple

import matplotlib.patches as patches
import matplotlib.pyplot as plt

from ..utils.styling import (
    ANNOTATION_TEXT_SIZE,
    TRACK_LABEL_COLOR,
    TRACK_LABEL_FONT_SIZE,
    get_distinct_colors,
)
from .base_track import BaseTrack

logger_ted_track = logging.getLogger(__name__)

PlottingOptionTED = Literal["full", "collapse"]


class TEDDomainsTrack(BaseTrack):
    """
    A track to display TED domain annotations.
    "full" mode: each original TED annotation (with its segments) is a row, colored uniquely, with CATH label.
    "collapse" mode: all segments from all annotations are merged and drawn in one row.
    """

    def __init__(
        self,
        ted_annotations: List[Dict[str, Any]],
        label: str = "TED Domains",
        plotting_option: PlottingOptionTED = "collapse",
        domain_height: float = 0.2,  # Height for individual segments/rows
        domain_spacing: float = 0.1,  # Vertical spacing in "full" mode
        default_color: str = "mediumpurple",  # For "collapse" mode or if palette runs out
        show_domain_labels: bool = True,  # For CATH labels in "full" mode
        domain_label_offset: float = 0.015,
        full_mode_padding: float = 0.15,
        collapse_mode_padding: float = 0.02,
    ):
        super().__init__(
            label=label,
            height=0,
            padding=full_mode_padding
            if plotting_option == "full"
            else collapse_mode_padding,
        )

        self.ted_annotations_orig = ted_annotations
        self.plotting_option = plotting_option
        self.domain_height = (
            domain_height  # Used as sub-track height in full, or bar height in collapse
        )
        self.domain_spacing = domain_spacing
        self.default_color = default_color
        self.show_domain_labels = show_domain_labels
        self.domain_label_offset = domain_label_offset
        self.full_mode_padding = full_mode_padding
        self.collapse_mode_padding = collapse_mode_padding

        # Colors for each original annotation entry in "full" mode
        self._annotation_colors = get_distinct_colors(
            len(self.ted_annotations_orig), len(self.ted_annotations_orig)
        )

        # Attributes for layout, populated by _layout_entries_for_view
        self._visible_annotations_full_mode: List[
            Dict[str, Any]
        ] = []  # Stores original annotation dicts that are visible
        self._parsed_segments_for_full_mode: Dict[
            int, List[Tuple[int, int]]
        ] = {}  # Key: original index of annotation
        self._merged_ted_regions_view: List[Dict[str, int]] = []

        if not self.ted_annotations_orig:
            self.height = self.domain_height  # Min height for "No data"

    def _parse_chopping_string(self, chopping_str: str) -> List[Tuple[int, int]]:
        segments = []
        if not chopping_str or not isinstance(chopping_str, str):
            return segments
        for part in chopping_str.split("_"):
            try:
                s_str, e_str = part.split("-")
                s, e = int(s_str), int(e_str)
                if s <= e:
                    segments.append((s, e))
                else:
                    logger_ted_track.warning(
                        f"Invalid segment in chopping '{part}' (start > end)."
                    )
            except ValueError:
                logger_ted_track.warning(f"Could not parse segment '{part}'.")
        return segments

    def _calculate_merged_ted_regions(
        self, all_segments: List[Tuple[int, int]]
    ) -> List[Dict[str, int]]:
        """Merges a flat list of (start, end) segments."""
        if not all_segments:
            return []

        # Sort by start, then by end
        sorted_segments = sorted(all_segments, key=lambda x: (x[0], x[1]))

        merged = []
        current_start, current_end = -1, -1

        for seg_start, seg_end in sorted_segments:
            if current_start == -1:
                current_start, current_end = seg_start, seg_end
            else:
                if seg_start <= current_end + 1:  # Overlap or adjacent
                    current_end = max(current_end, seg_end)
                else:  # New region
                    merged.append(
                        {"startIndex": current_start, "endIndex": current_end}
                    )
                    current_start, current_end = seg_start, seg_end
        if current_start != -1:
            merged.append({"startIndex": current_start, "endIndex": current_end})
        return merged

    def _layout_entries_for_view(self, view_start_aa: int, view_end_aa: int):
        if self.plotting_option == "collapse":
            self.padding = self.collapse_mode_padding
        else:
            self.padding = self.full_mode_padding

        self._visible_annotations_full_mode = []
        self._parsed_segments_for_full_mode = defaultdict(list)
        all_segments_for_collapse = []

        for i, annotation in enumerate(self.ted_annotations_orig):
            chopping_str = annotation.get("chopping")
            segments = self._parse_chopping_string(chopping_str)

            annotation_is_visible = False
            current_annotation_visible_segments = []
            for seg_start, seg_end in segments:
                # Check if segment overlaps with the view
                if not (seg_end < view_start_aa or seg_start > view_end_aa):
                    annotation_is_visible = True
                    current_annotation_visible_segments.append((seg_start, seg_end))
                    if (
                        self.plotting_option == "collapse"
                    ):  # Collect all segments for collapse mode
                        all_segments_for_collapse.append((seg_start, seg_end))

            if annotation_is_visible and self.plotting_option == "full":
                self._visible_annotations_full_mode.append(annotation)
                self._parsed_segments_for_full_mode[i] = (
                    current_annotation_visible_segments
                )

        if self.plotting_option == "full":
            num_lanes = len(self._visible_annotations_full_mode)
            if num_lanes > 0:
                self.height = (num_lanes * self.domain_height) + (
                    max(0, num_lanes - 1) * self.domain_spacing
                )
            else:
                self.height = self.domain_height  # For "No data" message

        elif self.plotting_option == "collapse":
            self._merged_ted_regions_view = self._calculate_merged_ted_regions(
                all_segments_for_collapse
            )
            self.height = self.domain_height  # Single row for collapsed view

        if (
            not self.ted_annotations_orig
            or (
                self.plotting_option == "full"
                and not self._visible_annotations_full_mode
            )
            or (
                self.plotting_option == "collapse" and not self._merged_ted_regions_view
            )
        ):
            self.height = self.domain_height

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

        no_data_to_draw = False
        if self.plotting_option == "full" and not self._visible_annotations_full_mode:
            no_data_to_draw = True
        if self.plotting_option == "collapse" and not self._merged_ted_regions_view:
            no_data_to_draw = True
        if not self.ted_annotations_orig:
            no_data_to_draw = True

        if no_data_to_draw:
            message_y_center = content_block_y_start + self.height / 2
            ax.text(
                (view_start_aa + view_end_aa) / 2,
                message_y_center,
                "No TED domains in view",
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
        # Iterate through original annotations that are determined to be visible
        for i, original_annotation_idx in enumerate(
            range(len(self.ted_annotations_orig))
        ):
            # Check if this original annotation is among the visible ones for full mode
            # This requires mapping index `i` here to the index in `self._visible_annotations_full_mode`
            # or more simply, iterate `self._visible_annotations_full_mode` and get their original index if needed for color.
            # Let's iterate based on the order of _visible_annotations_full_mode

            if i >= len(self._visible_annotations_full_mode):
                break  # Should not happen if logic is correct

            annotation_info = self._visible_annotations_full_mode[i]
            # Find original index for consistent coloring if _visible_annotations_full_mode doesn't store it
            # This is a bit complex. Simpler: color based on index in _visible_annotations_full_mode
            color = self._annotation_colors[
                self.ted_annotations_orig.index(annotation_info)
                % len(self._annotation_colors)
            ]

            lane_y_start = content_block_y_start + i * (
                self.domain_height + self.domain_spacing
            )

            # Get segments for this specific original annotation
            # The key for _parsed_segments_for_full_mode should be the original index
            original_idx_of_current_visible_annotation = (
                self.ted_annotations_orig.index(annotation_info)
            )
            segments_for_this_annotation = self._parsed_segments_for_full_mode.get(
                original_idx_of_current_visible_annotation, []
            )

            has_visible_segments_in_lane = False
            for seg_start, seg_end in segments_for_this_annotation:
                plot_s = max(seg_start, view_start_aa)
                plot_e = min(seg_end, view_end_aa)

                if plot_s < plot_e:
                    has_visible_segments_in_lane = True
                    bar_width = plot_e - plot_s
                    if bar_width <= 0 and plot_e == plot_s:
                        bar_width = 0.8
                    rect_x = (
                        plot_s - 0.4
                        if bar_width == 0.8 and (plot_e - plot_s) == 0
                        else plot_s
                    )

                    rect = patches.Rectangle(
                        (rect_x, lane_y_start),
                        bar_width,
                        self.domain_height,
                        linewidth=0.5,
                        edgecolor="black",
                        facecolor=color,
                        alpha=0.7,
                        clip_on=True,
                    )
                    ax.add_patch(rect)

            if self.show_domain_labels and has_visible_segments_in_lane:
                cath_label_raw = annotation_info.get("cath_label", "")
                display_label = (
                    "Not assigned" if cath_label_raw == "-" else cath_label_raw
                )

                if display_label:  # Only draw if there's something to display
                    label_y_center = lane_y_start + self.domain_height / 2
                    ax.text(
                        1.0 + self.domain_label_offset,
                        label_y_center,
                        display_label,
                        ha="left",
                        va="center",
                        fontsize=ANNOTATION_TEXT_SIZE,
                        color=color,
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
        # In collapse mode, all segments are drawn in a single row with default_color
        # The _merged_ted_regions_view already contains merged segments from *all* visible annotations

        lane_y_start = content_block_y_start  # Only one conceptual lane

        for region in self._merged_ted_regions_view:
            region_start = region["startIndex"]
            region_end = region["endIndex"]

            plot_s = max(region_start, view_start_aa)
            plot_e = min(region_end, view_end_aa)

            if plot_s < plot_e:
                bar_width = plot_e - plot_s
                if bar_width <= 0 and plot_e == plot_s:
                    bar_width = 0.8
                rect_x = (
                    plot_s - 0.4
                    if bar_width == 0.8 and (plot_e - plot_s) == 0
                    else plot_s
                )

                rect = patches.Rectangle(
                    (rect_x, lane_y_start),
                    bar_width,
                    self.domain_height,
                    linewidth=0.5,
                    edgecolor="black",
                    facecolor=self.default_color,
                    alpha=0.75,
                    clip_on=True,
                )
                ax.add_patch(rect)
