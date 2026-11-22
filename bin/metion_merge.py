#!/usr/bin/env python3
"""Merge per-comparison METION tables into a single report.

Writes:
  metion_summary.tsv   one row per region x modification, one column block per comparison
  metion_stats.tsv     run statistics for every comparison
  metion_report.xlsx   both tables plus the significant region sites of every comparison
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def merge_summaries(paths: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(p, sep="\t") for p in paths]
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame(columns=["region", "mod"])
    long = pd.concat(frames, ignore_index=True)
    wide = long.pivot_table(index=["region", "mod"], columns="comparison",
                            values=["n_hyper", "n_hypo"], fill_value=0, aggfunc="sum")
    wide.columns = [f"{comp}_{metric.replace('n_', '')}" for metric, comp in wide.columns]
    wide = wide[sorted(wide.columns)]
    return wide.reset_index()


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--summaries", nargs="+", type=Path, required=True)
    p.add_argument("--stats", nargs="+", type=Path, required=True)
    p.add_argument("--sites", nargs="*", type=Path, default=[],
                   help="*.region_sites.tsv files to add as sheets")
    args = p.parse_args(argv)

    summary = merge_summaries(args.summaries)
    stats = pd.concat([pd.read_csv(s, sep="\t") for s in args.stats], ignore_index=True)
    summary.to_csv("metion_summary.tsv", sep="\t", index=False)
    stats.to_csv("metion_stats.tsv", sep="\t", index=False)

    with pd.ExcelWriter("metion_report.xlsx") as xl:
        stats.to_excel(xl, sheet_name="stats", index=False)
        summary.to_excel(xl, sheet_name="regions_summary", index=False)
        for path in sorted(args.sites):
            name = path.name.replace(".region_sites.tsv", "")[:31]
            pd.read_csv(path, sep="\t").to_excel(xl, sheet_name=name, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
