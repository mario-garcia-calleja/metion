#!/usr/bin/env python3
"""Filter and annotate single-site differential methylation calls from `modkit dmr pair`.

Input is the BED produced by `modkit dmr pair` (single-site mode, 16 columns) and a
BED file of regions of interest (genes, promoters, operons...). The script:

1. parses the per-site calls and the modification code of each site;
2. applies coverage, effect-size and FDR (Benjamini-Hochberg) filters;
3. annotates significant sites with the regions they fall in, including a
   configurable flank (upstream only when the region strand is known);
4. writes genome-wide and region-level tables.

`--legacy` reproduces the original TFG behaviour (score > 0, gene body only,
no coverage/FDR filter) so results can be compared with the first version.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

LOG = logging.getLogger("metion_dmr")

DMR_COLUMNS = [
    "chrom", "start", "end", "name", "score", "strand",
    "a_counts", "a_total", "b_counts", "b_total",
    "a_mod_percentages", "b_mod_percentages",
    "a_pct_modified", "b_pct_modified",
    "map_pvalue", "effect_size",
]

MOD_NAMES = {"a": "6mA", "m": "5mC", "h": "5hmC", "21839": "4mC"}


# --------------------------------------------------------------------------- io
def read_dmr(path: Path) -> pd.DataFrame:
    """Read a `modkit dmr pair` single-site BED (with or without a '#' header)."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", names=DMR_COLUMNS,
                     dtype={"chrom": str, "a_counts": str, "b_counts": str})
    if df.empty:
        LOG.warning("%s contains no sites", path)
    # a_counts looks like "a:3" or, when several codes share a base, "m:2,21839:0"
    codes = df["a_counts"].str.replace(r":\d+", "", regex=True).str.split(",")
    df["mod_code"] = codes.str.join(",")
    df["mod"] = codes.apply(lambda cs: "+".join(MOD_NAMES.get(c, c) for c in cs))
    return df


def read_regions(path: Path) -> pd.DataFrame:
    """Read a BED3+ file; name (col 4) and strand (col 6) are optional."""
    raw = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
    if raw.shape[1] < 3:
        raise ValueError(f"{path}: regions BED needs at least 3 columns")
    regions = pd.DataFrame({
        "chrom": raw[0],
        "start": raw[1].astype(int),
        "end": raw[2].astype(int),
    })
    regions["region"] = raw[3] if raw.shape[1] > 3 else (
        regions["chrom"] + ":" + raw[1] + "-" + raw[2])
    regions["region_strand"] = raw[5] if raw.shape[1] > 5 else "."
    regions["region_strand"] = regions["region_strand"].fillna(".")
    return regions


# ----------------------------------------------------------------------- stats
def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaNs are kept as NaN)."""
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n == 0:
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


# ------------------------------------------------------------------ annotation
def annotate(sites: pd.DataFrame, regions: pd.DataFrame, flank: int,
             body_only: bool = False) -> pd.DataFrame:
    """Attach every overlapping region to each site (one row per site-region pair).

    A site overlaps a region if it falls in the region body, or (unless
    `body_only`) within `flank` bp upstream of it. When the region strand is
    unknown ('.') the flank is applied on both sides.
    """
    hits = []
    for reg in regions.itertuples(index=False):
        on_chrom = sites[sites["chrom"] == reg.chrom]
        if on_chrom.empty:
            continue
        if body_only or flank == 0:
            lo, hi = reg.start, reg.end
        elif reg.region_strand == "+":
            lo, hi = reg.start - flank, reg.end
        elif reg.region_strand == "-":
            lo, hi = reg.start, reg.end + flank
        else:
            lo, hi = reg.start - flank, reg.end + flank
        sel = on_chrom[(on_chrom["start"] >= lo) & (on_chrom["end"] <= hi)].copy()
        if sel.empty:
            continue
        in_body = (sel["start"] >= reg.start) & (sel["end"] <= reg.end)
        sel["region"] = reg.region
        sel["region_part"] = np.where(in_body, "body", "flank")
        hits.append(sel)
    if not hits:
        return sites.iloc[0:0].assign(region=pd.Series(dtype=str),
                                      region_part=pd.Series(dtype=str))
    return pd.concat(hits, ignore_index=True)


def summarise(annotated: pd.DataFrame) -> pd.DataFrame:
    """Count hyper/hypo-methylated sites per region and modification."""
    if annotated.empty:
        return pd.DataFrame(columns=["region", "mod", "n_sites", "n_hyper", "n_hypo"])
    grp = annotated.groupby(["region", "mod"])
    out = grp.size().rename("n_sites").to_frame()
    out["n_hyper"] = grp["direction"].apply(lambda s: (s == "hyper").sum())
    out["n_hypo"] = grp["direction"].apply(lambda s: (s == "hypo").sum())
    return out.reset_index().sort_values(["region", "mod"])


# ------------------------------------------------------------------------ main
def run(args: argparse.Namespace) -> None:
    sites = read_dmr(args.dmr)
    regions = read_regions(args.regions)
    LOG.info("%s: %d sites, %d regions", args.prefix, len(sites), len(regions))

    sites["fdr"] = bh_fdr(sites["map_pvalue"].to_numpy())
    # modkit reports percentages as 0-100; effect_size = a - b (fractions)
    sites["direction"] = np.where(sites["b_pct_modified"] > sites["a_pct_modified"],
                                  "hyper", "hypo")

    if args.legacy:
        sig = sites[sites["score"] > 0].copy()
        annotated = annotate(sig, regions, flank=0, body_only=True)
    else:
        cov_ok = (sites["a_total"] >= args.min_coverage) & (sites["b_total"] >= args.min_coverage)
        sig = sites[cov_ok
                    & (sites["fdr"] <= args.fdr)
                    & (sites["effect_size"].abs() >= args.min_effect)].copy()
        annotated = annotate(sig, regions, flank=args.flank)

    LOG.info("%s: %d significant sites, %d in regions of interest",
             args.prefix, len(sig), len(annotated))

    out_cols = ["chrom", "start", "end", "strand", "mod", "a_total", "b_total",
                "a_pct_modified", "b_pct_modified", "effect_size", "score",
                "map_pvalue", "fdr", "direction"]
    sig[out_cols].to_csv(f"{args.prefix}.significant_sites.tsv", sep="\t", index=False)
    annotated[["region", "region_part"] + out_cols].to_csv(
        f"{args.prefix}.region_sites.tsv", sep="\t", index=False)
    summary = summarise(annotated)
    summary.insert(0, "comparison", args.prefix)
    summary.to_csv(f"{args.prefix}.region_summary.tsv", sep="\t", index=False)

    stats = pd.DataFrame([{
        "comparison": args.prefix,
        "mode": "legacy" if args.legacy else "fdr",
        "sites_tested": len(sites),
        "sites_significant": len(sig),
        "sites_in_regions": annotated[["chrom", "start", "strand"]].drop_duplicates().shape[0],
        **{f"significant_{m}": int((sig["mod"] == m).sum()) for m in sorted(sites["mod"].unique())},
    }])
    stats.to_csv(f"{args.prefix}.stats.tsv", sep="\t", index=False)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--dmr", type=Path, required=True, help="modkit dmr pair output BED")
    p.add_argument("--regions", type=Path, required=True, help="BED of regions of interest")
    p.add_argument("--prefix", required=True, help="output prefix (comparison name)")
    p.add_argument("--fdr", type=float, default=0.05, help="max BH-adjusted p-value [0.05]")
    p.add_argument("--min-effect", type=float, default=0.2,
                   help="min |effect size| (difference in fraction modified) [0.2]")
    p.add_argument("--min-coverage", type=int, default=5,
                   help="min valid coverage in both samples [5]")
    p.add_argument("--flank", type=int, default=200,
                   help="bp upstream of each region also annotated (both sides if strand unknown) [200]")
    p.add_argument("--legacy", action="store_true",
                   help="reproduce the original TFG filter (score > 0, gene body only)")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))
    return 0


if __name__ == "__main__":
    sys.exit(main())
