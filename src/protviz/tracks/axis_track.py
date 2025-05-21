# --- src/protviz/tracks/axis_track.py ---
from typing import Optional

import matplotlib.pyplot as plt

from ..utils.styling import TRACK_LABEL_COLOR
from .base_track import BaseTrack


class AxisTrack(BaseTrack):
    """
    A simple track to draw the main axis representing the protein sequence.
    Adapts to zoomed views.
    """

    def __init__(
        self,
        sequence_length: int,  # This is the total sequence length
        label: str = "Sequence",
        height: float = 0.2,
        tick_interval: Optional[
            int
        ] = None,  # Optional: user can suggest, or auto-calculate
        major_tick_color: str = "black",
        minor_tick_color: str = "gray",
        text_color: str = "black",
        text_size: int = 8,
    ):
        super().__init__(label=label, height=height)
        self.total_sequence_length = sequence_length  # Store total length for context
        self.tick_interval_hint = tick_interval
        self.major_tick_color = major_tick_color
        self.minor_tick_color = (
            minor_tick_color  # Not implemented yet, but good for future
        )
        self.text_color = text_color
        self.text_size = text_size

    def _calculate_tick_interval(self, view_range_length: int) -> int:
        """Automatically determines a reasonable tick interval based on view range."""
        if self.tick_interval_hint:
            return self.tick_interval_hint

        if view_range_length <= 0:
            return 100  # Fallback

        # Aim for roughly 5-10 major ticks in the view
        if view_range_length <= 50:
            return 5
        if view_range_length <= 100:
            return 10
        if view_range_length <= 250:
            return 25
        if view_range_length <= 500:
            return 50
        if view_range_length <= 1000:
            return 100
        if view_range_length <= 2000:
            return 200
        if view_range_length <= 5000:
            return 500
        return 1000

    def draw(
        self,
        ax: plt.Axes,
        protein_sequence_length: int,
        track_y_bottom: float,
        view_start_aa: int,
        view_end_aa: int,
    ):
        content_y_center = track_y_bottom + self.padding + self.height / 2
        view_range_length = view_end_aa - view_start_aa + 1

        ax.plot(
            [view_start_aa, view_end_aa],
            [content_y_center, content_y_center],
            color=self.major_tick_color,
            lw=2,
            solid_capstyle="butt",
        )

        current_tick_interval = self._calculate_tick_interval(view_range_length)
        first_tick = (view_start_aa // current_tick_interval) * current_tick_interval
        if first_tick < view_start_aa:
            first_tick += current_tick_interval

        for i in range(first_tick, view_end_aa + 1, current_tick_interval):
            if i == 0 and view_start_aa > 1:
                continue
            if i < view_start_aa:
                continue

            tick_half_height = (
                self.height * 0.2
            )  # Reduced tick height: use a smaller fraction of self.height, e.g., height * 0.2
            ax.plot(
                [i, i],
                [
                    content_y_center - tick_half_height,
                    content_y_center + tick_half_height,
                ],
                color=self.major_tick_color,
                lw=1,
            )

            # Adjusted text position: closer to the axis line
            text_y_offset = self.height * 0.3  # Offset below the center line
            ax.text(
                i,
                content_y_center - text_y_offset,
                str(i),
                ha="center",
                va="top",
                fontsize=self.text_size,
                color=self.text_color,
            )
        if self.label:
            ax.text(
                0,
                content_y_center,
                self.label,
                ha="right",
                va="center",
                fontsize=self.text_size + 1,
                color=TRACK_LABEL_COLOR,
                transform=ax.get_yaxis_transform(),
                clip_on=False,
            )
