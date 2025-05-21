# --- src/protviz/utils/styling.py ---
# Purpose: To define common styling elements like color palettes and font settings.
from typing import List, Tuple

# Default font sizes
DEFAULT_FONT_SIZE = 8
ANNOTATION_TEXT_SIZE = 7
ROW_LABEL_TEXT_SIZE = 8
TRACK_LABEL_FONT_SIZE = 9

# Default colors
TRACK_LABEL_COLOR = "gray"
DEFAULT_FEATURE_COLOR = "royalblue"
ERROR_COLOR = "red"
WARNING_COLOR = "orange"


# pLDDT Score Colors (mimicking AlphaFold standard)
PLDDT_VERY_HIGH_COLOR = "#0052d6"  # Dark Blue (>90)
PLDDT_CONFIDENT_COLOR = "#65cbf3"  # Light Blue (70-90)
PLDDT_LOW_COLOR = "#FFDB13"  # Yellow (50-70)
PLDDT_VERY_LOW_COLOR = "#FF7D45"  # Orange (<50)

# AlphaMissense Pathogenicity Colors (example gradient: benign to pathogenic)
# This is a simple example; a more nuanced colormap might be better.
AM_BENIGN_COLOR = "#2166ac"  # Blue (e.g., < 0.34)
AM_AMBIGUOUS_COLOR = "#a8a9ac"  # Gray (e.g., 0.34 - 0.564)
AM_PATHOGENIC_COLOR = "#b2182b"  # Red (e.g., > 0.564)


# A simple categorical color palette (seaborn's "pastel" first few colors)
# (R, G, B) tuples, normalized to 0-1
_DEFAULT_PALETTE_RGB = [
    (0.6, 0.78, 1.0),  # Light Blue
    (1.0, 0.7, 0.6),  # Light Orange/Salmon
    (0.6, 1.0, 0.6),  # Light Green
    (0.9, 0.6, 0.8),  # Light Purple/Pink
    (1.0, 1.0, 0.6),  # Light Yellow
    (0.7, 0.9, 0.9),  # Light Cyan/Teal
]


def get_color_palette(num_colors: int) -> List[Tuple[float, float, float]]:
    """
    Returns a list of distinct colors.
    Cycles through the default palette if more colors are requested than available.
    """
    if num_colors <= 0:
        return []
    palette = []
    for i in range(num_colors):
        palette.append(_DEFAULT_PALETTE_RGB[i % len(_DEFAULT_PALETTE_RGB)])
    return palette


def get_plddt_color(score: float) -> str:
    if score >= 90:
        return PLDDT_VERY_HIGH_COLOR
    elif score >= 70:
        return PLDDT_CONFIDENT_COLOR
    elif score >= 50:
        return PLDDT_LOW_COLOR
    else:
        return PLDDT_VERY_LOW_COLOR


def get_alphamissense_color(avg_pathogenicity: float) -> str:
    # Example thresholds based on typical AlphaMissense interpretations
    # These thresholds might need adjustment based on desired visualization.
    if avg_pathogenicity < 0.34:  # Likely benign
        return AM_BENIGN_COLOR
    elif avg_pathogenicity <= 0.564:  # Ambiguous
        return AM_AMBIGUOUS_COLOR
    else:  # Likely pathogenic
        return AM_PATHOGENIC_COLOR


_BASE_DISTINCT_COLORS = [  # Store the base list
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#aec7e8",
    "#ffbb78",
    "#98df8a",
    "#ff9896",
    "#c5b0d5",
    "#c49c94",
    "#f7b6d2",
    "#c7c7c7",
    "#dbdb8d",
    "#9edae5",
]


def get_distinct_colors(num_colors: int, offset: int = 0) -> List[str]:
    """
    Provides a list of distinct hex color strings, with an optional offset.

    Args:
        num_colors (int): The number of distinct colors required.
        offset (int): An offset to start picking colors from the base list.
                      The offset will be taken modulo the length of the base color list.

    Returns:
        List[str]: A list of hex color strings.
    """
    if num_colors <= 0:
        return []

    base_len = len(_BASE_DISTINCT_COLORS)
    actual_offset = offset % base_len if base_len > 0 else 0

    palette = []
    for i in range(num_colors):
        palette.append(_BASE_DISTINCT_COLORS[(actual_offset + i) % base_len])
    return palette
