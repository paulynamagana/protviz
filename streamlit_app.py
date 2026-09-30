"""Protviz web app.

Run locally with:  streamlit run streamlit_app.py
Deployed on Streamlit Community Cloud, which picks up this file automatically.
"""

import io
import re

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from matplotlib.colors import is_color_like

from protviz import plot_protein_tracks
from protviz.data_retrieval import (
    AFDBClient,
    InterProClient,
    PDBeClient,
    TEDClient,
    get_protein_sequence_length,
)
from protviz.tracks import (
    AlphaFoldTrack,
    AxisTrack,
    CustomTrack,
    InterProTrack,
    LigandInteractionTrack,
    PDBTrack,
    TEDDomainsTrack,
)

UNIPROT_ID_PATTERN = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)

st.set_page_config(page_title="Protviz", page_icon="🧬", layout="wide")


# --- Data fetching (cached so re-plotting or zooming doesn't re-download) ---


@st.cache_resource
def get_clients():
    return {
        "pdbe": PDBeClient(),
        "afdb": AFDBClient(),
        "ted": TEDClient(),
        "interpro": InterProClient(),
    }


@st.cache_data(ttl=86400, show_spinner=False)
def fetch_length(uniprot_id):
    return get_protein_sequence_length(uniprot_id)


@st.cache_data(ttl=86400, show_spinner=False)
def fetch_pdb_coverage(uniprot_id):
    return get_clients()["pdbe"].get_pdb_coverage(uniprot_id)


@st.cache_data(ttl=86400, show_spinner=False)
def fetch_ligands(uniprot_id):
    return get_clients()["pdbe"].get_pdb_ligand_interactions(uniprot_id)


@st.cache_data(ttl=86400, show_spinner=False)
def fetch_alphafold(uniprot_id):
    return get_clients()["afdb"].get_alphafold_data(
        uniprot_id, requested_data_types=["plddt", "alphamissense"]
    )


@st.cache_data(ttl=86400, show_spinner=False)
def fetch_ted(uniprot_id):
    return get_clients()["ted"].get_TED_annotations(uniprot_id)


@st.cache_data(ttl=86400, show_spinner=False)
def fetch_pfam(uniprot_id):
    return get_clients()["interpro"].get_pfam_annotations(uniprot_id)


@st.cache_data(ttl=86400, show_spinner=False)
def fetch_cath(uniprot_id):
    return get_clients()["interpro"].get_cathgene3d_annotations(uniprot_id)


def safe_fetch(fetcher, uniprot_id, source_name):
    """Fetch data for one track; on failure warn and return None so other tracks still plot."""
    try:
        return fetcher(uniprot_id)
    except Exception:
        st.warning(
            f"Could not load {source_name} data right now, so that track was skipped."
        )
        return None


def custom_rows_to_annotations(df):
    """Turn the editable table into CustomTrack annotation dicts, skipping incomplete rows."""

    def text(value, default):
        return (
            default
            if value is None or pd.isna(value) or not str(value).strip()
            else str(value).strip()
        )

    annotations = []
    for row in df.to_dict("records"):
        start = row.get("Start")
        if start is None or pd.isna(start):
            continue
        end = row.get("End")
        entry = {
            "label": text(row.get("Label"), ""),
            "row_label": text(row.get("Row"), "Custom"),
            "color": text(row.get("Color"), "royalblue"),
        }
        if not is_color_like(entry["color"]):
            st.warning(
                f"'{entry['color']}' isn't a color I recognise, so blue was used instead."
            )
            entry["color"] = "royalblue"
        if end is None or pd.isna(end) or int(end) == int(start):
            entry["position"] = int(start)
            entry["display_type"] = "marker"
        else:
            entry["start"], entry["end"] = sorted((int(start), int(end)))
        annotations.append(entry)
    return annotations


# --- Sidebar: inputs ---

with st.sidebar:
    st.title("🧬 Protviz")
    st.caption("Plot protein annotations along the sequence.")

    uniprot_id = (
        st.text_input(
            "UniProt ID",
            value="O15245",
            help="For example P00533 (EGFR) or O15245. Find IDs at uniprot.org.",
        )
        .strip()
        .upper()
    )

    st.subheader("Tracks")
    show_pdb = st.checkbox("PDB structure coverage", value=True)
    pdb_mode = st.radio(
        "PDB display",
        ["collapse", "full"],
        horizontal=True,
        disabled=not show_pdb,
        help="'collapse' shows overall coverage; 'full' shows each PDB entry.",
    )
    show_ligands = st.checkbox("Ligand binding sites", value=False)
    show_plddt = st.checkbox("AlphaFold confidence (pLDDT)", value=True)
    show_am = st.checkbox(
        "AlphaMissense pathogenicity",
        value=False,
        help="Only available for human proteins.",
    )
    show_ted = st.checkbox("TED domains", value=False)
    show_pfam = st.checkbox("Pfam domains", value=False)
    show_cath = st.checkbox("CATH-Gene3D domains", value=False)

    st.subheader("Plot size")
    figure_width = st.slider("Width", 6, 24, 14)

# --- Main area ---

st.header("Protein annotation viewer")

if not uniprot_id:
    st.info("Enter a UniProt ID in the sidebar to get started.")
    st.stop()

if not UNIPROT_ID_PATTERN.match(uniprot_id):
    st.error(
        f"'{uniprot_id}' doesn't look like a UniProt ID. "
        "IDs look like P00533 or O15245."
    )
    st.stop()

with st.spinner(f"Looking up {uniprot_id} in UniProt..."):
    try:
        seq_length = fetch_length(uniprot_id)
    except (ValueError, KeyError):
        st.error(
            f"UniProt has no active entry with a sequence for '{uniprot_id}'. "
            "Check the ID and try again."
        )
        st.stop()
    except Exception:
        st.error("UniProt could not be reached right now. Please try again shortly.")
        st.stop()

view_start, view_end = st.slider(
    "Zoom to region (amino acids)",
    min_value=1,
    max_value=seq_length,
    value=(1, seq_length),
)

with st.expander("Add your own annotations"):
    st.caption(
        "One row per feature. Leave **End** empty to mark a single residue. "
        "Rows with the same **Row** name are drawn on the same line."
    )
    custom_df = st.data_editor(
        pd.DataFrame(
            {
                "Start": pd.Series([], dtype="Int64"),
                "End": pd.Series([], dtype="Int64"),
                "Label": pd.Series([], dtype="str"),
                "Row": pd.Series([], dtype="str"),
                "Color": pd.Series([], dtype="str"),
            }
        ),
        num_rows="dynamic",
        width="stretch",
        column_config={
            "Start": st.column_config.NumberColumn(min_value=1, max_value=seq_length),
            "End": st.column_config.NumberColumn(min_value=1, max_value=seq_length),
            "Color": st.column_config.TextColumn(
                help="A color name (red, teal) or hex code (#ff8800). Defaults to blue."
            ),
        },
        key=f"custom_{uniprot_id}",
    )

tracks = [AxisTrack(sequence_length=seq_length, label="Sequence")]

with st.spinner("Fetching annotations..."):
    if show_pdb:
        data = safe_fetch(fetch_pdb_coverage, uniprot_id, "PDB")
        if data is not None:
            tracks.append(
                PDBTrack(pdb_data=data, label="PDB", plotting_option=pdb_mode)
            )

    if show_ligands:
        data = safe_fetch(fetch_ligands, uniprot_id, "ligand")
        if data is not None:
            tracks.append(
                LigandInteractionTrack(
                    interaction_data=data,
                    label="Ligands",
                    plotting_option="collapse",
                    show_ligand_labels=True,
                )
            )

    af_options = [
        opt for opt, on in (("plddt", show_plddt), ("alphamissense", show_am)) if on
    ]
    if af_options:
        data = safe_fetch(fetch_alphafold, uniprot_id, "AlphaFold")
        if data is not None:
            if show_am and not data.get("alphamissense"):
                st.info(
                    f"No AlphaMissense data for {uniprot_id} (it covers human proteins only)."
                )
            tracks.append(
                AlphaFoldTrack(
                    afdb_data=data,
                    plotting_options=af_options,
                    main_label="",
                    plddt_label="pLDDT",
                    alphamissense_label="AlphaMissense",
                    sub_track_height=0.1,
                    sub_track_spacing=0.05,
                )
            )

    if show_ted:
        data = safe_fetch(fetch_ted, uniprot_id, "TED")
        if data is not None:
            tracks.append(
                TEDDomainsTrack(
                    ted_annotations=data, label="TED", plotting_option="full"
                )
            )

    if show_pfam:
        data = safe_fetch(fetch_pfam, uniprot_id, "Pfam")
        if data is not None:
            tracks.append(
                InterProTrack(
                    domain_data=data, database_name_for_label="Pfam", label="Pfam"
                )
            )

    if show_cath:
        data = safe_fetch(fetch_cath, uniprot_id, "CATH-Gene3D")
        if data is not None:
            tracks.append(
                InterProTrack(
                    domain_data=data,
                    database_name_for_label="CATH-Gene3D",
                    label="CATH-Gene3D",
                )
            )

    custom_annotations = custom_rows_to_annotations(custom_df)
    if custom_annotations:
        tracks.append(
            CustomTrack(
                annotation_data=custom_annotations,
                show_row_labels=True,
                show_ann_labels=True,
                ann_height=0.1,
                ann_spacing=0.05,
            )
        )

if len(tracks) == 1:
    st.info("Tick at least one track in the sidebar to see annotations.")

try:
    fig = plot_protein_tracks(
        protein_id=uniprot_id,
        sequence_length=seq_length,
        tracks=tracks,
        figure_width=figure_width,
        view_start_aa=view_start,
        view_end_aa=view_end,
        show=False,
    )
except Exception:
    st.error(
        "Something went wrong while drawing the plot. "
        "Try turning off the most recently added track or custom annotation."
    )
    st.stop()

st.pyplot(fig, width="stretch")

buffer = io.BytesIO()
fig.savefig(buffer, format="png", dpi=300, bbox_inches="tight")
plt.close(fig)

st.download_button(
    "Download PNG (300 dpi)",
    data=buffer.getvalue(),
    file_name=f"{uniprot_id}_protviz.png",
    mime="image/png",
)
