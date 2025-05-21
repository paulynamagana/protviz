# src/protviz/tracks/interpro_track.py
import logging
from typing import Any, Dict, List, Literal, Optional, Tuple

import matplotlib.patches as patches
import matplotlib.pyplot as plt

from ..utils.styling import (
    ANNOTATION_TEXT_SIZE,
    DEFAULT_FEATURE_COLOR,
    TRACK_LABEL_COLOR,
    TRACK_LABEL_FONT_SIZE,
    get_distinct_colors,
)
from .base_track import BaseTrack

logger_ip_track = logging.getLogger(__name__)

PlottingOptionInterPro = Literal["full", "collapse"]

# Mapping InterPro entry types to single characters for labels
INTERPRO_TYPE_CHAR_MAP = {
    "domain": "D",
    "family": "F",
    "homologous_superfamily": "H",
    "repeat": "R",
    "site": "S",
    "unknown": "U",
    "ptm": "P",
}


class InterProTrack(BaseTrack):
    """
    A track to display domain annotations from InterPro member databases (e.g., Pfam, CATH-Gene3D).
    Supports "full" mode (each member DB signature in its own row with detailed labels)
    and "collapse" mode (all segments merged into a single row).
    """

    def __init__(
        self,
        domain_data: List[Dict[str, Any]],
        database_name_for_label: str,
        label: str = "Domains",
        plotting_option: PlottingOptionInterPro = "full",
        domain_height: Optional[float] = 0.2,
        domain_spacing: Optional[float] = 0.05,
        default_color: str = DEFAULT_FEATURE_COLOR,  # For "collapse" mode
        show_domain_labels: bool = True,  # For right-side labels in "full" mode
        domain_label_offset: float = 0.015,
        label_line_spacing: float = 0.008,  # Vertical spacing between lines of right-side label
        full_mode_padding: Optional[float] = 0.05,
        collapse_mode_padding: Optional[float] = 0.05,
    ):
        super().__init__(
            label=label,
            height=0,
            padding=full_mode_padding
            if plotting_option == "full"
            else collapse_mode_padding,
        )

        self.domain_data_orig = domain_data
        self.database_name_for_label = (
            database_name_for_label.upper()
        )  # e.g. PFAM, CATHGENE3D
        self.plotting_option = plotting_option
        self.domain_height = domain_height
        self.domain_spacing = domain_spacing
        self.default_color = default_color
        self.show_domain_labels = show_domain_labels
        self.domain_label_offset = domain_label_offset
        self.label_line_spacing = label_line_spacing  # For multi-line right labels
        self.full_mode_padding = full_mode_padding
        self.collapse_mode_padding = collapse_mode_padding

        # For "full" mode layout
        self._visible_domains_full_mode: List[Dict[str, Any]] = []
        self._domain_colors_full_mode: Dict[str, str] = {}  # Keyed by domain accession

        # For "collapse" mode layout
        self._merged_domain_regions_view: List[Dict[str, int]] = []

        if not self.domain_data_orig:
            self.height = self.domain_height

    def _parse_domain_locations(
        self, annotation: Dict[str, Any]
    ) -> List[Tuple[int, int]]:
        """Extracts and flattens all (start, end) segments for a given domain annotation."""
        all_segments = []
        for loc in annotation.get("locations", []):
            if (
                isinstance(loc, dict)
                and loc.get("start") is not None
                and loc.get("end") is not None
            ):
                try:
                    all_segments.append((int(loc["start"]), int(loc["end"])))
                except ValueError:
                    logger_ip_track.warning(
                        f"Invalid location in domain {annotation.get('accession')}: {loc}"
                    )
        return all_segments

    def _calculate_merged_regions(
        self, all_segments_flat: List[Tuple[int, int]]
    ) -> List[Dict[str, int]]:
        if not all_segments_flat:
            return []
        sorted_segments = sorted(all_segments_flat, key=lambda x: (x[0], x[1]))
        merged = []
        current_start, current_end = -1, -1
        for seg_start, seg_end in sorted_segments:
            if current_start == -1:
                current_start, current_end = seg_start, seg_end
            else:
                if seg_start <= current_end + 1:
                    current_end = max(current_end, seg_end)
                else:
                    merged.append({"start": current_start, "end": current_end})
                    current_start, current_end = seg_start, seg_end
        if current_start != -1:
            merged.append({"start": current_start, "end": current_end})
        return merged

    def _layout_entries_for_view(self, view_start_aa: int, view_end_aa: int):
        if self.plotting_option == "collapse":
            self.padding = self.collapse_mode_padding
        else:
            self.padding = self.full_mode_padding

        self._visible_domains_full_mode = []
        all_segments_for_collapse = []

        # Filter original annotations for visibility
        temp_visible_domains = []
        for ann in self.domain_data_orig:
            segments = self._parse_domain_locations(ann)
            is_visible_in_view = any(
                not (seg_end < view_start_aa or seg_start > view_end_aa)
                for seg_start, seg_end in segments
            )
            if is_visible_in_view:
                temp_visible_domains.append(
                    ann
                )  # Store the whole annotation if any part is visible
                if self.plotting_option == "collapse":
                    for (
                        seg_start,
                        seg_end,
                    ) in segments:  # Collect all segments for merging
                        if not (
                            seg_end < view_start_aa or seg_start > view_end_aa
                        ):  # Double check visibility for segment
                            all_segments_for_collapse.append((seg_start, seg_end))

        self._visible_domains_full_mode = sorted(
            temp_visible_domains, key=lambda x: x.get("accession", "")
        )

        if self.plotting_option == "full":
            num_lanes = len(self._visible_domains_full_mode)
            self._domain_colors_full_mode = {
                ann.get("accession"): color
                for ann, color in zip(
                    self._visible_domains_full_mode,
                    get_distinct_colors(num_lanes, num_lanes),
                )
            }
            if num_lanes > 0:
                self.height = (num_lanes * self.domain_height) + (
                    max(0, num_lanes - 1) * self.domain_spacing
                )
            else:
                self.height = self.domain_height

        elif self.plotting_option == "collapse":
            self._merged_domain_regions_view = self._calculate_merged_regions(
                all_segments_for_collapse
            )
            self.height = self.domain_height

        if (
            not self.domain_data_orig
            or (self.plotting_option == "full" and not self._visible_domains_full_mode)
            or (
                self.plotting_option == "collapse"
                and not self._merged_domain_regions_view
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

        no_data = (
            (self.plotting_option == "full" and not self._visible_domains_full_mode)
            or (
                self.plotting_option == "collapse"
                and not self._merged_domain_regions_view
            )
            or not self.domain_data_orig
        )

        if no_data:
            msg_y = content_block_y_start + self.height / 2
            ax.text(
                (view_start_aa + view_end_aa) / 2,
                msg_y,
                f"No {self.database_name_for_label} data in view",
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
        for i, domain_annotation in enumerate(self._visible_domains_full_mode):
            domain_acc = domain_annotation.get("accession", "N/A")
            domain_name = domain_annotation.get("name", "")
            interpro_type = domain_annotation.get(
                "entry_type", "unknown"
            )  # Get InterPro type
            type_char = INTERPRO_TYPE_CHAR_MAP.get(interpro_type, "U")

            color = self._domain_colors_full_mode.get(domain_acc, self.default_color)
            lane_y_start = content_block_y_start + i * (
                self.domain_height + self.domain_spacing
            )

            segments = self._parse_domain_locations(domain_annotation)
            has_visible_segments = False
            for seg_start, seg_end in segments:
                plot_s, plot_e = (
                    max(seg_start, view_start_aa),
                    min(seg_end, view_end_aa),
                )
                if plot_s < plot_e:
                    has_visible_segments = True
                    bar_w = plot_e - plot_s
                    rect = patches.Rectangle(
                        (plot_s, lane_y_start),
                        bar_w,
                        self.domain_height,
                        linewidth=0.5,
                        edgecolor="black",
                        facecolor=color,
                        alpha=0.75,
                        clip_on=True,
                    )
                    ax.add_patch(rect)

            if self.show_domain_labels and has_visible_segments:
                label_y_center = lane_y_start + self.domain_height / 2
                # Line 1: Type char and Name
                label_text_line1 = f"{type_char}: {domain_name}"
                # Line 2: DB Name and Accession
                label_text_line2 = f"{self.database_name_for_label}: {domain_acc}"

                # Calculate y positions for two lines of text
                text_y1 = (
                    label_y_center + self.label_line_spacing * self.domain_height
                )  # Upper line
                text_y2 = (
                    label_y_center - self.label_line_spacing * self.domain_height
                )  # Lower line

                ax.text(
                    1.0 + self.domain_label_offset,
                    text_y1,
                    label_text_line1,
                    ha="left",
                    va="bottom",
                    fontsize=ANNOTATION_TEXT_SIZE - 1,
                    color=color,  # Smaller for two lines
                    transform=ax.get_yaxis_transform(),
                    clip_on=False,
                )
                ax.text(
                    1.0 + self.domain_label_offset,
                    text_y2,
                    label_text_line2,
                    ha="left",
                    va="top",
                    fontsize=ANNOTATION_TEXT_SIZE - 1,
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
        lane_y_start = content_block_y_start
        for region in self._merged_domain_regions_view:
            plot_s, plot_e = (
                max(region["start"], view_start_aa),
                min(region["end"], view_end_aa),
            )
            if plot_s < plot_e:
                bar_w = plot_e - plot_s
                rect = patches.Rectangle(
                    (plot_s, lane_y_start),
                    bar_w,
                    self.domain_height,
                    linewidth=0.5,
                    edgecolor="black",
                    facecolor=self.default_color,
                    alpha=0.8,
                    clip_on=True,
                )
                ax.add_patch(rect)
