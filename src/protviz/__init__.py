# src/protviz/__init__.py
from . import tracks  # Makes the tracks module available
from .core_plotting import plot_protein_tracks
from .data_retrieval import get_protein_sequence_length
from .tracks import AxisTrack, BaseTrack, PDBTrack

__version__ = "0.0.1"

__all__ = [
    "plot_protein_tracks",
    "get_protein_sequence_length",
    "tracks",
    "AxisTrack",
    "PDBTrack",
    "BaseTrack",
]
