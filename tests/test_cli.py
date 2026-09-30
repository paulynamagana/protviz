import pytest

from protviz.cli import (
    DEFAULT_TRACKS,
    TRACK_ORDER,
    UserError,
    build_parser,
    main,
    parse_region,
    resolve_tracks,
)


class TestParseRegion:
    @pytest.mark.parametrize(
        "text, expected",
        [
            ("700-1000", (700, 1000)),
            (" 700 - 1000 ", (700, 1000)),
            ("700:1000", (700, 1000)),
            ("1-2", (1, 2)),
        ],
    )
    def test_valid_regions(self, text, expected):
        assert parse_region(text) == expected

    @pytest.mark.parametrize("text", ["banana", "700", "", "abc-def", "700-"])
    def test_unparseable_regions_raise_user_error(self, text):
        with pytest.raises(UserError):
            parse_region(text)

    def test_start_must_be_positive(self):
        with pytest.raises(UserError, match="must be 1 or more"):
            parse_region("0-100")

    def test_reversed_region_suggests_the_swap(self):
        with pytest.raises(UserError, match=r"Try --region 100-900"):
            parse_region("900-100")


class TestResolveTracks:
    def test_no_selection_gives_the_defaults(self):
        assert resolve_tracks(None) == DEFAULT_TRACKS
        assert resolve_tracks([]) == DEFAULT_TRACKS

    def test_all_gives_every_track(self):
        assert resolve_tracks(["all"]) == TRACK_ORDER

    def test_accepts_commas_and_spaces(self):
        assert resolve_tracks(["pdb,pfam"]) == resolve_tracks(["pdb", "pfam"])

    def test_result_is_in_display_order_not_argument_order(self):
        assert resolve_tracks(["plddt", "pdb"]) == ["pdb", "plddt"]

    def test_names_are_case_insensitive(self):
        assert resolve_tracks(["PDB", "Pfam"]) == ["pdb", "pfam"]

    def test_duplicates_are_collapsed(self):
        assert resolve_tracks(["pdb", "pdb"]) == ["pdb"]

    def test_unknown_track_names_are_listed_back_to_the_user(self):
        with pytest.raises(UserError) as excinfo:
            resolve_tracks(["pdb", "nonsense"])
        message = str(excinfo.value)
        assert "nonsense" in message
        assert "--list-tracks" in message


class TestCommandLine:
    def test_list_tracks_exits_cleanly_and_names_every_track(self, capsys):
        assert main(["--list-tracks"]) == 0
        output = capsys.readouterr().out
        for name in TRACK_ORDER:
            assert name in output

    def test_missing_accession_is_an_error_with_guidance(self, capsys):
        assert main([]) == 1
        assert "uniprot.org" in capsys.readouterr().err

    def test_bad_region_fails_before_any_network_call(self, capsys):
        assert main(["P04637", "--region", "banana"]) == 1
        assert "START-END" in capsys.readouterr().err

    def test_defaults(self):
        args = build_parser().parse_args(["P04637"])
        assert args.uniprot_id == "P04637"
        assert args.width == 12.0
        assert args.dpi == 300
        assert not args.show
        assert not args.detail
