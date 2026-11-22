#!/usr/bin/env python3
"""Generate the small synthetic dataset used by `-profile test` and the unit tests.

Everything is simulated (no real sequencing data is distributed with the repository):
  reference.fasta  20 kb random genome ('chrT') containing GATC motifs
  regions.bed      three 'genes' on chrT
  *.bed.gz         bedMethyl files for 4 samples (control/treated x 2), bgzipped + tabix

Run from the repository root:  python tests/make_test_data.py
Requires bgzip and tabix on PATH.
"""
from __future__ import annotations

import random
import subprocess
from pathlib import Path

OUT = Path(__file__).parent / "data"
CHROM = "chrT"
SEED = 42


def make_reference(rng: random.Random, length: int = 20_000) -> str:
    seq = list(rng.choices("ACGT", k=length))
    for pos in range(100, length - 100, 250):  # plant a GATC every 250 bp
        seq[pos:pos + 4] = list("GATC")
    return "".join(seq)


def write_fasta(path: Path, name: str, seq: str) -> None:
    lines = [seq[i:i + 70] for i in range(0, len(seq), 70)]
    path.write_text(f">{name}\n" + "\n".join(lines) + "\n")


def bedmethyl_row(pos: int, strand: str, code: str, cov: int, n_mod: int) -> str:
    pct = 100 * n_mod / cov
    return "\t".join(map(str, [
        CHROM, pos, pos + 1, code, cov, strand, pos, pos + 1, "255,0,0",
        cov, f"{pct:.2f}", n_mod, cov - n_mod, 0, 0, 0, 0, 0]))


def make_bedmethyl(rng: random.Random, ref: str, treated: bool, path: Path) -> None:
    rows = []
    for pos, base in enumerate(ref):
        # the same sites are sequenced in every sample: A/C bases at every 5th position
        # plus the adenine of every GATC motif
        if base not in "AC" or not (pos % 5 == 0 or ref[pos - 1:pos + 3] == "GATC"):
            continue
        code = "a" if base == "A" else "m"
        cov = rng.randint(8, 20)
        in_gatc = ref[pos - 1:pos + 3] == "GATC"
        frac = 0.95 if in_gatc else 0.03
        # treated samples lose Dam methylation inside gene_B (positions 9000-11000)
        if treated and in_gatc and 9_000 <= pos < 11_000:
            frac = 0.05
        n_mod = sum(rng.random() < frac for _ in range(cov))
        rows.append(bedmethyl_row(pos, "+", code, cov, n_mod))
    path.write_text("\n".join(rows) + "\n")
    subprocess.run(["bgzip", "-f", str(path)], check=True)
    subprocess.run(["tabix", "-f", "-p", "bed", f"{path}.gz"], check=True)


def main() -> None:
    rng = random.Random(SEED)
    OUT.mkdir(parents=True, exist_ok=True)

    ref = make_reference(rng)
    write_fasta(OUT / "reference.fasta", CHROM, ref)
    (OUT / "regions.bed").write_text(
        f"{CHROM}\t2000\t4000\tgene_A\t0\t+\n"
        f"{CHROM}\t9000\t11000\tgene_B\t0\t-\n"
        f"{CHROM}\t15000\t16500\tgene_C\t0\t.\n")

    # METION: two control/treated pairs starting from bedMethyl
    sheet = ["sample,barcode,bam,bedmethyl"]
    for sample, treated in [("ctrl_1", False), ("trt_1", True), ("ctrl_2", False), ("trt_2", True)]:
        make_bedmethyl(rng, ref, treated, OUT / f"{sample}.bed")
        sheet.append(f"{sample},,,{sample}.bed.gz")
    (OUT / "samplesheet.csv").write_text("\n".join(sheet) + "\n")
    (OUT / "comparisons.csv").write_text(
        "comparison,control,treated\npair1,ctrl_1,trt_1\npair2,ctrl_2,trt_2\n")


if __name__ == "__main__":
    main()
