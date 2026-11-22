"""Unit tests for the METION analysis scripts."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def load(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[path.stem] = module
    spec.loader.exec_module(module)
    return module


dmr = load(ROOT / "bin" / "metion_dmr.py")


# ------------------------------------------------------------------- METION
def test_bh_fdr_matches_reference_values():
    p = np.array([0.01, 0.04, 0.03, 0.2, np.nan])
    adj = dmr.bh_fdr(p)
    # R: p.adjust(c(0.01, 0.04, 0.03, 0.2), "BH") -> 0.04 0.0533 0.0533 0.2
    np.testing.assert_allclose(adj[:4], [0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.2])
    assert np.isnan(adj[4])


def _sites(positions):
    return pd.DataFrame({"chrom": "c", "start": positions, "end": [p + 1 for p in positions]})


def test_annotate_uses_upstream_flank_by_strand():
    regions = pd.DataFrame({"chrom": ["c", "c"], "start": [1000, 5000], "end": [2000, 6000],
                            "region": ["plus", "minus"], "region_strand": ["+", "-"]})
    sites = _sites([900, 1500, 2100, 4900, 6100])
    out = dmr.annotate(sites, regions, flank=200)
    got = set(zip(out["region"], out["start"], out["region_part"]))
    assert got == {("plus", 900, "flank"), ("plus", 1500, "body"), ("minus", 6100, "flank")}


def test_annotate_unknown_strand_flanks_both_sides():
    regions = pd.DataFrame({"chrom": ["c"], "start": [1000], "end": [2000],
                            "region": ["g"], "region_strand": ["."]})
    out = dmr.annotate(_sites([850, 2150, 2300]), regions, flank=200)
    assert sorted(out["start"]) == [850, 2150]


def _dmr_file(tmp_path: Path) -> Path:
    rows = [
        # strong, well-covered change inside the region -> significant
        ["c", 1500, 1501, ".", 9.0, "+", "a:2", 20, "a:19", 20, "a:10.00", "a:95.00", 10.0, 95.0, 1e-6, 0.85],
        # large change but low coverage -> filtered
        ["c", 1600, 1601, ".", 2.0, "+", "a:0", 2, "a:2", 2, "a:0.00", "a:100.00", 0.0, 100.0, 1e-3, 1.0],
        # no change -> filtered (but score > 0 keeps it in legacy mode)
        ["c", 1700, 1701, ".", 0.1, "+", "m:1", 20, "m:1", 20, "m:5.00", "m:5.00", 5.0, 5.0, 0.9, 0.0],
    ]
    path = tmp_path / "x.dmr.bed"
    pd.DataFrame(rows).to_csv(path, sep="\t", header=False, index=False)
    return path


def test_dmr_filters_and_legacy_mode(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    regions = tmp_path / "r.bed"
    regions.write_text("c\t1000\t2000\tgeneX\t0\t+\n")
    dmr_bed = _dmr_file(tmp_path)

    dmr.main(["--dmr", str(dmr_bed), "--regions", str(regions), "--prefix", "new"])
    sig = pd.read_csv("new.significant_sites.tsv", sep="\t")
    assert list(sig["start"]) == [1500]
    assert sig.loc[0, "mod"] == "6mA" and sig.loc[0, "direction"] == "hyper"
    summary = pd.read_csv("new.region_summary.tsv", sep="\t")
    assert summary.to_dict("records") == [
        {"comparison": "new", "region": "geneX", "mod": "6mA", "n_sites": 1, "n_hyper": 1, "n_hypo": 0}]

    dmr.main(["--dmr", str(dmr_bed), "--regions", str(regions), "--prefix", "old", "--legacy"])
    assert len(pd.read_csv("old.significant_sites.tsv", sep="\t")) == 3


def test_read_dmr_names_combined_cytosine_codes(tmp_path):
    row = ["c", 10, 11, ".", 1.0, "+", "m:2,21839:1", 10, "m:5,21839:3", 10,
           "m:20.00,21839:10.00", "m:50.00,21839:30.00", 0.3, 0.8, 0.01, -0.5]
    path = tmp_path / "d.bed"
    pd.DataFrame([row]).to_csv(path, sep="\t", header=False, index=False)
    assert dmr.read_dmr(path).loc[0, "mod"] == "5mC+4mC"
