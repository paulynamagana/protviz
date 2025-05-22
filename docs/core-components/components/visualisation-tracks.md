---
layout: default
title: Visualisation Tracks
parent: Core Components
nav_order: 2
description: "Understand the Visualisation Tracks."
---

# Visualisation Tracks
<button class="btn js-toggle-dark-mode">Dark mode</button>

<script>
const toggleDarkMode = document.querySelector('.js-toggle-dark-mode');

jtd.addEvent(toggleDarkMode, 'click', function(){
  if (jtd.getTheme() === 'dark') {
    jtd.setTheme('light');
    toggleDarkMode.textContent = 'Dark mode';
  } else {
    jtd.setTheme('dark');
    toggleDarkMode.textContent = 'Light mode';
  }
});
</script>
---
Once you have retrieved your data, **Visualisation Tracks** are used to display it on the protein plot. Each track is a layer showing a specific type of annotation. Protviz offers a variety of specialised tracks, and you can combine multiple tracks to create a comprehensive view of your protein.
{: .fs-6 .fw-300 }

All tracks inherit from a `BaseTrack` class, which provides common functionalities like handling labels and height. When you create a plot using `plot_protein_tracks()`, you provide a list of these track objects. They are then drawn in order, typically from the bottom of the plot upwards.

Here's a detailed look at the available track types:

1.  **`AxisTrack`**
    * **Purpose**: Displays the main horizontal axis representing the protein sequence. It includes numbered tick marks and an optional label, providing essential context for all other annotations.
    * **Key Parameters**:
        * `sequence_length` (int): The total length of the protein. This is crucial for scaling the axis correctly.
        * `label` (str, optional): A descriptive label for the axis (e.g., "Sequence", "Protein Length"). Defaults to "Sequence".
        * `tick_interval` (int, optional): Allows you to suggest an interval for the major tick marks. If omitted, Protviz calculates a sensible default based on the viewable sequence range.
    * **Usage Example**:
        ```python
        from protviz.tracks import AxisTrack

        # Assuming seq_len is already fetched
        # seq_len = get_protein_sequence_length("P0DTC2")
        axis = AxisTrack(sequence_length=seq_len, label="SARS-CoV-2 Spike Protein")
        ```

2.  **`PDBTrack`**
    * **Purpose**: Visualises regions of the protein that are covered by structural entries in the Protein Data Bank (PDB).
    * **Data Input**: Expects a list of dictionaries, typically from `PDBeClient().get_pdb_coverage()`, where each dictionary details a PDB entry's start and end on the UniProt sequence.
    * **Key Parameters**:
        * `pdb_data` (List[Dict]): The PDB coverage data.
        * `label` (str, optional): Label for the track (e.g., "PDB Coverage").
        * `plotting_option` (str, optional):
            * `"full"`: Each PDB entry is drawn in a separate lane, attempting to label them with their PDB ID if space allows. This helps distinguish individual structural contributions.
            * `"collapse"` (default): Merges all overlapping PDB covered regions into a single, continuous bar, giving an overview of overall structural coverage.
        * `bar_height` (float, optional): The height of the PDB bars.
        * `color` (str, optional): The fill colour for the PDB bars.
    * **Usage Example**:
        ```python
        from protviz.tracks import PDBTrack
        # Assuming pdbe_client is initialised and pdb_coverage_data is fetched
        # pdb_coverage_data = pdbe_client.get_pdb_coverage("P00533")
        pdb_track_full = PDBTrack(pdb_data=pdb_coverage_data, label="PDB Entries (Full)", plotting_option="full")
        pdb_track_summary = PDBTrack(pdb_data=pdb_coverage_data, label="PDB Coverage (Summary)", plotting_option="collapse", color="lightcoral")
        ```

3.  **`LigandInteractionTrack`**
    * **Purpose**: Highlights sites on the protein that are known to interact with ligands, based on data from PDB entries.
    * **Data Input**: Expects a list of dictionaries from `PDBeClient().get_pdb_ligand_interactions()`, where each dictionary represents a ligand and its interacting residue ranges on the UniProt sequence.
    * **Key Parameters**:
        * `interaction_data` (List[Dict]): The ligand interaction data.
        * `label` (str, optional): Label for the track (e.g., "Ligand Binding Sites").
        * `plotting_option` (str, optional):
            * `"full"` (default): Each unique ligand is shown in its own lane, with its binding sites represented as bars. Ligand IDs are typically displayed.
            * `"collapse"`: All ligand binding sites from all ligands are merged into a single representative bar.
        * `show_ligand_labels` (bool, optional): If `True` (default), displays ligand IDs next to their respective lanes in "full" mode.
        * `site_height` (float, optional): Height of the bars representing binding sites.
        * `default_color` (str, optional): Colour used in "collapse" mode, or if distinct colours for ligands run out.
    * **Usage Example**:
        ```python
        from protviz.tracks import LigandInteractionTrack
        # Assuming pdbe_client is initialised and ligand_data is fetched
        # ligand_data = pdbe_client.get_pdb_ligand_interactions("P00533")
        ligand_track = LigandInteractionTrack(
            interaction_data=ligand_data,
            label="Ligand Interactions",
            plotting_option="full",
            show_ligand_labels=True
        )
        ```

4.  **`TEDDomainsTrack`**
    * **Purpose**: Displays structural domain annotations sourced from the TED database. TED annotations often include CATH domain classifications and defined segment boundaries ("chopping").
    * **Data Input**: Expects a list of dictionaries, usually from `TEDClient().get_TED_annotations()`.
    * **Key Parameters**:
        * `ted_annotations` (List[Dict]): The TED annotation data.
        * `label` (str, optional): Label for the track (e.g., "TED Domains").
        * `plotting_option` (str, optional):
            * `"full"`: Each original TED annotation (which can consist of multiple segments) is drawn in its own lane, coloured uniquely, and can be labelled with its CATH classification.
            * `"collapse"` (default): All domain segments from all TED annotations are merged and drawn as a single bar.
        * `show_domain_labels` (bool, optional): If `True` (default), displays CATH labels for each lane in "full" mode.
        * `domain_height` (float, optional): Height of the domain bars/segments.
    * **Usage Example**:
        ```python
        from protviz.tracks import TEDDomainsTrack
        # Assuming ted_client is initialised and ted_domain_data is fetched
        # ted_domain_data = ted_client.get_TED_annotations("O15245")
        ted_track_detailed = TEDDomainsTrack(
            ted_annotations=ted_domain_data,
            label="TED Domains (Detailed View)",
            plotting_option="full",
            show_domain_labels=True
        )
        ```

5.  **`AlphaFoldTrack`**
    * **Purpose**: Visualises prediction metrics from the AlphaFold Database, primarily the per-residue pLDDT confidence scores and (average) AlphaMissense pathogenicity scores.
    * **Data Input**: Expects a dictionary from `AFDBClient().get_alphafold_data()`, where keys are data types (e.g., `"plddt"`, `"alphamissense"`) and values are lists of corresponding score entries.
    * **Key Parameters**:
        * `afdb_data` (Dict[str, List[Dict]]): The AlphaFold data.
        * `plotting_options` (List[str], optional): A list specifying which data types to plot from the `afdb_data` dictionary (e.g., `["plddt"]`, `["alphamissense"]`, or `["plddt", "alphamissense"]`). If `None`, it attempts to plot all supported types found in the data.
        * `main_label` (str, optional): An overall label for the group of AlphaFold sub-tracks.
        * `plddt_label` (str, optional): Custom label for the pLDDT sub-track. Defaults to "pLDDT".
        * `alphamissense_label` (str, optional): Custom label for the AlphaMissense sub-track. Defaults to "AlphaMissense (avg)".
        * `sub_track_height` (float, optional): Height for each individual row (pLDDT, AlphaMissense).
    * **Usage Example**:
        ```python
        from protviz.tracks import AlphaFoldTrack
        # Assuming afdb_client is initialised and alphafold_prediction_data is fetched
        # alphafold_prediction_data = afdb_client.get_alphafold_data("P0DTC2", requested_data_types=["plddt", "alphamissense"])
        alphafold_viewer_track = AlphaFoldTrack(
            afdb_data=alphafold_prediction_data,
            plotting_options=["plddt", "alphamissense"],
            main_label="AlphaFold Predictions",
            plddt_label="Model Confidence (pLDDT)",
            alphamissense_label="Pathogenicity (AM)"
        )
        ```

6.  **`InterProTrack`**
    * **Purpose**: Displays domain and feature annotations sourced from InterPro member databases like Pfam or CATH-Gene3D.
    * **Data Input**: Expects a list of domain annotation dictionaries, typically obtained from methods like `InterProClient().get_pfam_annotations()` or `InterProClient().get_cathgene3d_annotations()`.
    * **Key Parameters**:
        * `domain_data` (List[Dict]): The list of domain annotations.
        * `database_name_for_label` (str): The name of the member database (e.g., "Pfam", "CATH-Gene3D"). This is used in labels if `show_domain_labels` is active.
        * `label` (str, optional): A main label for the track.
        * `plotting_option` (str, optional):
            * `"full"` (default): Each unique domain signature from the member database is displayed in its own row. Labels on the right can show the domain's type (from InterPro), name, and its specific accession from the member database.
            * `"collapse"`: All domain segments are merged into a single representative bar.
        * `show_domain_labels` (bool, optional): If `True` (default), displays detailed labels on the right side in "full" mode.
        * `domain_height` (float, optional): The height of the domain bars.
    * **Usage Example**:
        ```python
        from protviz.tracks import InterProTrack
        # Assuming interpro_client is initialised and pfam_data is fetched
        # pfam_data = interpro_client.get_pfam_annotations("P04637")
        pfam_domains_track = InterProTrack(
            domain_data=pfam_data,
            database_name_for_label="Pfam",
            label="Pfam Domains",
            plotting_option="full"
        )
        ```

7.  **`CustomTrack`**
    * **Purpose**: Offers maximum flexibility by allowing you to plot your own, arbitrary annotations. This is ideal for data not covered by the specialised tracks, such as post-translational modification (PTM) sites, specific motifs, or user-defined regions of interest.
    * **Data Input**: A list of dictionaries, where each dictionary defines a custom annotation.
    * **Key Parameters for each annotation item in `annotation_data`**:
        * `start` (int) and `end` (int): For range-based annotations (displayed as bars by default).
        * OR `position` (int): For point-based annotations (displayed as markers by default). Protviz treats this as `start` = `end` = `position`.
        * `label` (str, optional): A specific label for this individual annotation, typically shown on the right.
        * `row_label` (str, optional): A label to group annotations into horizontal lanes, shown on the left. If multiple annotations share the same `row_label`, they will effectively be grouped, though each still gets its own visual lane in the current implementation.
        * `color` (str, optional): Colour for the annotation.
        * `display_type` (str, optional): Can be `"bar"` or `"marker"`. If not specified, it defaults to `"marker"` for point annotations (`start` == `end`) and `"bar"` for range annotations.
        * `marker_symbol` (str, optional): Matplotlib marker symbol (e.g., `"o"`, `"P"`, `"*"`), used if `display_type="marker"`.
        * `marker_size` (int, optional): Size of the marker.
    * **Key Parameters for the `CustomTrack` itself**:
        * `annotation_data` (List[Dict]): The list of custom annotation dictionaries.
        * `label` (str, optional): An overall label for the entire custom track group.
        * `show_row_labels` (bool, optional): Whether to display the `row_label`s on the left. Default is `True`.
        * `show_ann_labels` (bool, optional): Whether to display the individual annotation `label`s on the right. Default is `True`.
        * `ann_height` (float, optional): Height of each annotation's visual lane.
    * **Usage Example**:
        ```python
        from protviz.tracks import CustomTrack

        my_annotations = [
            {"position": 45, "label": "S45-P", "row_label": "Phosphorylation", "color": "blue", "display_type": "marker", "marker_symbol": "P"},
            {"start": 100, "end": 150, "label": "Binding Motif", "row_label": "Motifs", "color": "green"},
            {"position": 188, "label": "Y188-OH", "row_label": "Hydroxylation", "color": "purple", "display_type": "marker", "marker_symbol": "*"},
            {"start": 200, "end": 210, "label": "Signal Peptide", "row_label": "Processing", "color": "orange"}
        ]
        user_defined_track = CustomTrack(annotation_data=my_annotations, label="My Custom Features")
        ```

By combining these diverse tracks, you can construct highly informative and tailored visualisations of your protein data. The final appearance can be further tuned using parameters within `plot_protein_tracks()` such as figure width and zoom levels.