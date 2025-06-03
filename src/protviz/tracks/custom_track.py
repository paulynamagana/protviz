# --- src/protviz/tracks/custom_track.py
import logging
from typing import Any, Dict, List, Literal, Optional

import matplotlib.patches as patches
import matplotlib.pyplot as plt

from ..utils.styling import (
    ANNOTATION_TEXT_SIZE,
    DEFAULT_FEATURE_COLOR,
    TRACK_LABEL_COLOR,
    TRACK_LABEL_FONT_SIZE,
)
from .base_track import BaseTrack

logger_custom_track = logging.getLogger(__name__)

DisplayType = Literal["bar", "marker"]


class CustomTrack(BaseTrack):
    """
    A flexible track to display various types of custom annotations.
    Each annotation item gets its own horizontal lane.
    'row_label' can be used for left-side categorization labels.
    'label' is used for right-side specific annotation labels.
    """

    def __init__(
        self,
        annotation_data: List[Dict[str, Any]],
        label: str = "",  # Main label for the entire track group
        show_row_labels: bool = True,
        default_row_label: str = "",  # Default if 'row_label' is missing, effectively no row label
        ann_height: float = 0.1,  # Height of each annotation's lane
        ann_spacing: float = 0.05,  # Vertical spacing between annotation lanes
        default_color: str = DEFAULT_FEATURE_COLOR,
        show_ann_labels: bool = True,  # Show individual annotation labels (on the right)
        ann_label_offset: float = 0.015,  # Offset for right-side labels
        ann_label_size: int = ANNOTATION_TEXT_SIZE,
        row_label_size: int = ANNOTATION_TEXT_SIZE,
        default_display_type_range: DisplayType = "bar",
        default_display_type_point: Optional[DisplayType] = "marker",
        default_marker_symbol: str = "o",
        default_marker_size: Optional[int] = 5,
        padding: float = 0.1,
    ):
        super().__init__(label=label, height=0, padding=padding)

        self.annotation_data_orig = annotation_data
        self.show_row_labels = show_row_labels
        self.default_row_label = (
            default_row_label  # Used if an item doesn't have 'row_label'
        )
        self.ann_height = ann_height
        self.ann_spacing = ann_spacing
        self.default_color = default_color
        self.show_ann_labels = show_ann_labels
        self.ann_label_offset = ann_label_offset
        self.ann_label_size = ann_label_size
        self.row_label_size = row_label_size
        self.default_display_type_range = default_display_type_range
        self.default_display_type_point = default_display_type_point
        self.default_marker_symbol = default_marker_symbol
        self.default_marker_size = default_marker_size

        self.processed_annotations = self._process_custom_ann_data(
            self.annotation_data_orig
        )

        # Will store the list of annotation dicts that are currently visible
        self._visible_annotations_in_view: List[Dict[str, Any]] = []

        if not self.processed_annotations:
            self.height = self.ann_height  # Min height for "No data" message

    def _process_custom_ann_data(
        self, raw_data: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        processed_data = []
        if not isinstance(raw_data, list):
            logger_custom_track.error("Annotation data must be a list of dictionaries.")
            return []

        for i, item_orig in enumerate(raw_data):
            if not isinstance(item_orig, dict):
                logger_custom_track.warning(
                    f"Annotation item at index {i} is not a dictionary. Skipping."
                )
                continue

            item = item_orig.copy()

            if "position" in item:
                if "start" in item or "end" in item:
                    logger_custom_track.warning(
                        f"Annotation '{item.get('label', i)}' has both 'position' and 'start'/'end'. Using 'position'."
                    )
                try:
                    pos = int(item["position"])
                    item["start"] = pos
                    item["end"] = pos
                except ValueError:
                    logger_custom_track.warning(
                        f"Invalid 'position' for annotation '{item.get('label', i)}'. Skipping."
                    )
                    continue
            elif "start" not in item or "end" not in item:
                logger_custom_track.warning(
                    f"Annotation '{item.get('label', i)}' must have 'start'/'end' or 'position'. Skipping."
                )
                continue

            try:
                item["start"] = int(item["start"])
                item["end"] = int(item["end"])
            except ValueError:
                logger_custom_track.warning(
                    f"Invalid 'start'/'end' for annotation '{item.get('label', i)}'. Skipping."
                )
                continue

            if item["start"] > item["end"]:
                logger_custom_track.debug(
                    f"Swapping start/end for annotation '{item.get('label', i)}'."
                )
                item["start"], item["end"] = item["end"], item["start"]

            item.setdefault(
                "row_label", self.default_row_label
            )  # Default row label if not provided
            item.setdefault("color", self.default_color)
            item.setdefault(
                "label", ""
            )  # Specific label for this annotation item (for right side)

            is_point_annotation = item["start"] == item["end"]
            default_display = (
                self.default_display_type_point
                if is_point_annotation
                else self.default_display_type_range
            )
            item.setdefault("display_type", default_display)

            if item["display_type"] == "marker":
                item.setdefault("marker_symbol", self.default_marker_symbol)
                item.setdefault("marker_size", self.default_marker_size)

            processed_data.append(item)
        return processed_data

    def _layout_entries_for_view(self, view_start_aa: int, view_end_aa: int):
        """Filters annotations for visibility and calculates track height."""
        self._visible_annotations_in_view = []
        for ann in self.processed_annotations:
            s, e = ann["start"], ann["end"]
            if not (e < view_start_aa or s > view_end_aa):  # Check for any overlap
                self._visible_annotations_in_view.append(ann)

        # Height is based on the number of individual visible annotations, each getting a row
        num_visible_annotations = len(self._visible_annotations_in_view)
        if num_visible_annotations > 0:
            self.height = (num_visible_annotations * self.ann_height) + (
                max(0, num_visible_annotations - 1) * self.ann_spacing
            )
        else:
            self.height = self.ann_height  # For "No data" message

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

        if not self._visible_annotations_in_view:
            message_y_center = content_block_y_start + self.height / 2
            ax.text(
                (view_start_aa + view_end_aa) / 2,
                message_y_center,
                "No custom data in view",
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

        # Draw annotations from bottom up. Each annotation gets its own lane.
        for i, ann in enumerate(self._visible_annotations_in_view):
            # Y position for the bottom of the current annotation's lane
            lane_y_start = content_block_y_start + i * (
                self.ann_height + self.ann_spacing
            )
            lane_center_y = lane_y_start + self.ann_height / 2

            # Draw row label on the left (if applicable)
            row_label_text = ann.get("row_label", self.default_row_label)
            if (
                self.show_row_labels and row_label_text != self.default_row_label
            ):  # Only draw if not the default empty one
                # To avoid re-drawing for consecutive same labels, check against previous
                # For now, draw if present.
                ax.text(
                    -0.01,
                    lane_center_y,
                    row_label_text,  # Small x-offset from y-axis
                    ha="right",
                    va="center",
                    fontsize=self.row_label_size,
                    color=TRACK_LABEL_COLOR,
                    transform=ax.get_yaxis_transform(),
                    clip_on=False,
                )

            # Draw the annotation feature itself
            orig_start, orig_end = ann["start"], ann["end"]
            plot_s = max(orig_start, view_start_aa)
            plot_e = min(orig_end, view_end_aa)

            if plot_s > plot_e:
                continue  # Should be caught by _layout, but double check

            color = ann["color"]
            display_type = ann["display_type"]

            is_point = orig_start == orig_end  # Check if it's a point annotation
            bar_width = plot_e - plot_s  # Width of the bar for range annotations

            # Adjust bar width and x position for point annotations
            if is_point:
                bar_width = 0.8
                current_plot_x = plot_s - bar_width / 2
            else:
                bar_width = max(0.8, plot_e - plot_s)
                current_plot_x = plot_s

            if display_type == "bar":
                rect = patches.Rectangle(
                    (current_plot_x, lane_y_start),
                    bar_width,
                    self.ann_height,
                    facecolor=color,
                    edgecolor="black",
                    linewidth=0.5,
                    alpha=0.7,
                    clip_on=True,
                )
                ax.add_patch(rect)
            elif display_type == "marker":
                marker_pos = (plot_s + plot_e) / 2
                ax.plot(
                    marker_pos,
                    lane_center_y,
                    marker=ann.get("marker_symbol", self.default_marker_symbol),
                    color=color,
                    markersize=ann.get("marker_size", self.default_marker_size),
                    linestyle="None",
                    clip_on=True,
                )

            # Draw specific annotation label on the right side (unless it was a text_label type)
            specific_ann_label = ann.get("label", "")
            if (
                self.show_ann_labels
                and specific_ann_label
                and display_type != "text_label"
            ):
                ax.text(
                    1.0 + self.ann_label_offset,
                    lane_center_y,
                    specific_ann_label,
                    ha="left",
                    va="center",
                    fontsize=self.ann_label_size,
                    color=color,  # Use annotation color for its label
                    transform=ax.get_yaxis_transform(),
                    clip_on=False,
                )

        # Draw the main track label (for the whole group)
        if self.label:
            main_label_y_center = content_block_y_start + self.height / 2
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
