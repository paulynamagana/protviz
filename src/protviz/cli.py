"""Command-line interface for protviz.

Lets people make a protein annotation figure without writing any Python:

    protviz P00533
"""

import argparse
import logging
import sys
from typing import List, Optional, Tuple

# Order in which tracks are stacked in the figure, regardless of the order
# the user listed them on the command line.
TRACK_ORDER = ["pdb", "ligands", "ted", "pfam", "cath", "plddt", "alphamissense"]

TRACK_HELP = {
    "pdb": "Experimental structures from the PDB, shown as sequence coverage",
    "ligands": "Positions where a ligand binds the protein, from the PDB",
    "ted": "Structural domains predicted by TED",
    "pfam": "Protein families and domains from Pfam",
    "cath": "Structural domains from CATH-Gene3D",
    "plddt": "AlphaFold per-residue confidence (pLDDT)",
    "alphamissense": "AlphaMissense average pathogenicity (slow: downloads a large file)",
}

DEFAULT_TRACKS = ["pdb", "ligands", "ted", "pfam", "plddt"]


class UserError(Exception):
    """An error worth showing to the user as a plain sentence, with no traceback."""


def not_found_message(uniprot_id: str) -> str:
    return (
        f"UniProt has no entry for '{uniprot_id}'.\n\n"
        "protviz needs a UniProt accession, which looks like P00533 or Q9Y6K9. "
        "Gene names\nsuch as 'EGFR' or 'TP53' will not work.\n\n"
        "To find the accession, search for your protein at https://www.uniprot.org "
        "and copy\nthe short code shown in the 'Entry' column."
    )


def parse_region(text: str) -> Tuple[int, int]:
    """Turn "700-1000" into (700, 1000)."""
    separator = "-" if "-" in text else ":" if ":" in text else None
    if separator is None:
        raise UserError(
            f"Could not understand the region '{text}'. "
            "Write it as START-END, for example: --region 700-1000"
        )
    start_text, _, end_text = text.partition(separator)
    try:
        start, end = int(start_text.strip()), int(end_text.strip())
    except ValueError:
        raise UserError(
            f"Could not understand the region '{text}'. "
            "Both numbers must be whole numbers, for example: --region 700-1000"
        )
    if start < 1:
        raise UserError(
            f"The region start must be 1 or more, but you gave {start}. "
            "Amino acids are numbered starting at 1."
        )
    if start >= end:
        raise UserError(
            f"The region start ({start}) must be smaller than the end ({end}). "
            f"Try --region {end}-{start} instead."
        )
    return start, end


def resolve_tracks(requested: Optional[List[str]]) -> List[str]:
    """Validate the requested track names and put them in display order."""
    if not requested:
        return list(DEFAULT_TRACKS)

    names = [name.strip().lower() for item in requested for name in item.split(",")]
    names = [name for name in names if name]

    if "all" in names:
        return list(TRACK_ORDER)

    unknown = [name for name in names if name not in TRACK_HELP]
    if unknown:
        raise UserError(
            "Unknown track name{0}: {1}.\nAvailable tracks are: {2}.\n"
            "Run 'protviz --list-tracks' to see what each one shows.".format(
                "s" if len(unknown) > 1 else "",
                ", ".join(repr(name) for name in unknown),
                ", ".join(TRACK_ORDER),
            )
        )
    return [name for name in TRACK_ORDER if name in names]


def print_track_list() -> None:
    print("Tracks you can ask for with --tracks:\n")
    width = max(len(name) for name in TRACK_ORDER)
    for name in TRACK_ORDER:
        default_marker = "  (on by default)" if name in DEFAULT_TRACKS else ""
        print(f"  {name:<{width}}  {TRACK_HELP[name]}{default_marker}")
    print("\nUse --tracks all to include every track.")


def fetch_tracks(uniprot_id, seq_length, wanted, detail):
    """Fetch data for each requested track and build the track objects.

    A database that is down or has nothing for this protein produces a warning
    and a skipped track, never a crash.
    """
    from .data_retrieval import AFDBClient, InterProClient, PDBeClient, TEDClient
    from .tracks import (
        AlphaFoldTrack,
        AxisTrack,
        InterProTrack,
        LigandInteractionTrack,
        PDBTrack,
        TEDDomainsTrack,
    )

    mode = "full" if detail else "collapse"
    tracks = [AxisTrack(sequence_length=seq_length, label="Sequence")]

    def attempt(description, build):
        """Run one track builder, reporting progress and surviving failures."""
        print(f"  Fetching {description}...", file=sys.stderr)
        try:
            track = build()
        except Exception as error:
            print(
                f"  Skipped {description}: {error}",
                file=sys.stderr,
            )
            return
        if track is not None:
            tracks.append(track)

    if "pdb" in wanted:

        def build_pdb():
            data = PDBeClient().get_pdb_coverage(uniprot_id)
            if not data:
                print("  No PDB structures found for this protein.", file=sys.stderr)
                return None
            return PDBTrack(
                pdb_data=data,
                label="PDB structures",
                plotting_option=mode,
                color="skyblue",
            )

        attempt("PDB structure coverage", build_pdb)

    if "ligands" in wanted:

        def build_ligands():
            data = PDBeClient().get_pdb_ligand_interactions(uniprot_id)
            if not data:
                print("  No ligand binding sites found.", file=sys.stderr)
                return None
            return LigandInteractionTrack(
                interaction_data=data,
                label="Ligand binding",
                plotting_option=mode,
                show_ligand_labels=detail,
            )

        attempt("ligand binding sites", build_ligands)

    if "ted" in wanted:

        def build_ted():
            data = TEDClient().get_TED_annotations(uniprot_id)
            if not data:
                print("  No TED domains found.", file=sys.stderr)
                return None
            return TEDDomainsTrack(
                ted_annotations=data,
                label="TED domains",
                plotting_option=mode,
            )

        attempt("TED domains", build_ted)

    for name, label, fetch_name in [
        ("pfam", "Pfam domains", "get_pfam_annotations"),
        ("cath", "CATH-Gene3D domains", "get_cathgene3d_annotations"),
    ]:
        if name not in wanted:
            continue

        def build_interpro(label=label, fetch_name=fetch_name, name=name):
            data = getattr(InterProClient(), fetch_name)(uniprot_id)
            if not data:
                print(f"  No {label} found.", file=sys.stderr)
                return None
            return InterProTrack(
                domain_data=data,
                database_name_for_label=name.upper(),
                label=label,
                plotting_option="full" if detail else "collapse",
            )

        attempt(label, build_interpro)

    alphafold_options = [name for name in ("plddt", "alphamissense") if name in wanted]
    if alphafold_options:

        def build_alphafold():
            data = AFDBClient().get_alphafold_data(
                uniprot_id, requested_data_types=alphafold_options
            )
            if not any(data.get(option) for option in alphafold_options):
                print("  No AlphaFold data found for this protein.", file=sys.stderr)
                return None
            return AlphaFoldTrack(
                afdb_data=data,
                plotting_options=alphafold_options,
                main_label="AlphaFold",
            )

        description = " and ".join(alphafold_options)
        if "alphamissense" in alphafold_options:
            description += " (this one can take a minute)"
        attempt(f"AlphaFold {description}", build_alphafold)

    return tracks


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="protviz",
        description=(
            "Draw what is known about a protein - structures, domains, ligand "
            "binding sites and AlphaFold confidence - as a single figure."
        ),
        epilog=(
            "Examples:\n"
            "  protviz P00533                        Save P00533.png with the usual tracks\n"
            "  protviz P00533 -o egfr.png            Choose the file name\n"
            "  protviz P00533 --tracks pdb,pfam      Pick just the tracks you want\n"
            "  protviz P00533 --region 700-1000      Zoom into amino acids 700 to 1000\n"
            "  protviz P00533 --show                 Open a window instead of saving\n"
            "  protviz --list-tracks                 See every available track\n"
            "\nA UniProt accession looks like P00533 or Q9Y6K9. You can find one by\n"
            "searching your protein's name at https://www.uniprot.org"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "uniprot_id",
        nargs="?",
        help="UniProt accession of the protein to draw, for example P00533",
    )
    parser.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help="Where to save the image. Defaults to <UNIPROT_ID>.png in this folder.",
    )
    parser.add_argument(
        "-t",
        "--tracks",
        nargs="+",
        metavar="NAME",
        help=(
            "Which tracks to draw, separated by spaces or commas. "
            f"Defaults to: {', '.join(DEFAULT_TRACKS)}. Use 'all' for everything."
        ),
    )
    parser.add_argument(
        "-r",
        "--region",
        metavar="START-END",
        help="Zoom into part of the protein, for example 700-1000.",
    )
    parser.add_argument(
        "--detail",
        action="store_true",
        help="Show every entry on its own row instead of merging them into one bar.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open the figure in a window instead of saving it to a file.",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=12.0,
        metavar="INCHES",
        help="Width of the figure in inches (default: 12).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        metavar="N",
        help="Resolution of the saved image (default: 300).",
    )
    parser.add_argument(
        "--list-tracks",
        action="store_true",
        help="List every available track and what it shows, then exit.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show detailed progress messages, useful when something goes wrong.",
    )
    return parser


def run(args) -> int:
    if args.list_tracks:
        print_track_list()
        return 0

    if not args.uniprot_id:
        raise UserError(
            "Please give a UniProt accession, for example:\n"
            "  protviz P00533\n\n"
            "You can find the accession for your protein by searching its name at "
            "https://www.uniprot.org"
        )

    wanted = resolve_tracks(args.tracks)
    region = parse_region(args.region) if args.region else None

    saving = bool(args.out) or not args.show
    if saving:
        import matplotlib

        matplotlib.use("Agg")

    from .core_plotting import plot_protein_tracks
    from .data_retrieval import get_protein_sequence_length

    uniprot_id = args.uniprot_id.strip().upper()

    import requests

    print(f"Looking up {uniprot_id} on UniProt...", file=sys.stderr)
    try:
        seq_length = get_protein_sequence_length(uniprot_id)
    except (ValueError, KeyError):
        raise UserError(not_found_message(uniprot_id))
    except requests.exceptions.HTTPError as error:
        # UniProt answers 400 for text that is not accession-shaped, 404 for one
        # that is well formed but unused.
        status = getattr(error.response, "status_code", None)
        if status in (400, 404):
            raise UserError(not_found_message(uniprot_id))
        raise UserError(
            f"UniProt returned an error (HTTP {status}) for '{uniprot_id}'.\n"
            "This is usually temporary - please try again in a moment."
        )
    except Exception as error:
        raise UserError(
            f"Could not reach UniProt to look up '{uniprot_id}'.\n"
            "Check your internet connection and try again.\n\n"
            f"Technical detail: {error}"
        )

    print(f"{uniprot_id} is {seq_length} amino acids long.", file=sys.stderr)

    if region and region[0] > seq_length:
        raise UserError(
            f"You asked for region {region[0]}-{region[1]}, but {uniprot_id} is only "
            f"{seq_length} amino acids long."
        )

    tracks = fetch_tracks(uniprot_id, seq_length, wanted, args.detail)

    if len(tracks) <= 1:
        raise UserError(
            f"None of the requested tracks had any data for {uniprot_id}, so there is "
            "nothing to draw.\nTry a different protein, or run with --tracks all to "
            "widen the search."
        )

    output_path = args.out or (f"{uniprot_id}.png" if saving else None)

    plot_protein_tracks(
        protein_id=uniprot_id,
        sequence_length=seq_length,
        tracks=tracks,
        figure_width=args.width,
        view_start_aa=region[0] if region else None,
        view_end_aa=region[1] if region else None,
        save_path=output_path,
        show=args.show,
        dpi=args.dpi,
    )
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s - %(name)s - %(message)s",
    )

    try:
        return run(args)
    except UserError as error:
        print(f"\nError: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nStopped.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
