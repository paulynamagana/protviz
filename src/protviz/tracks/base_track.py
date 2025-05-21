# --- src/protviz/tracks/base_track.py ---
from abc import ABC, abstractmethod

import matplotlib.pyplot as plt


class BaseTrack(ABC):
    def __init__(self, label: str = None, height: float = 0.5, padding: float = 0.1):
        self.label = label
        self.height = height
        self.padding = padding

    @abstractmethod
    def draw(
        self,
        ax: plt.Axes,
        protein_sequence_length: int,
        track_y_bottom: float,
        view_start_aa: int,
        view_end_aa: int,
    ):
        pass

    def get_total_height(self) -> float:
        """Returns the total height this track will occupy including top/bottom padding."""
        return self.height + self.padding * 1.5
