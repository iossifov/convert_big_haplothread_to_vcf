#!/usr/bin/env python

import pandas as pd
import sys
import gzip

in_data_dir = "/mnt/wigclust8/data/safe/autism/HT/AGRE_WG38_AGRE859"
coll = "AGRE859"
chrom = "chr18"
chunk_size = 10000
out_dir = "./chunks"

assert len(sys.argv) in [1, 6]
if len(sys.argv) == 6:
    in_data_dir, coll, chrom, chunk_size, out_dir = sys.argv[1:]
    chunk_size = int(chunk_size)

print(f"Working on {in_data_dir=}, {coll=}, {chrom=}, {chunk_size=}, {out_dir=}...")

in_hpth_file = f"{in_data_dir}/{coll}-{chrom}-hpth.txt.gz"
in_pos_file = f"{in_data_dir}/{coll}-{chrom}-pos.txt"
meta_data_file = f"{in_data_dir}/metaData.txt"

PS = pd.read_table(in_pos_file, sep='\t', header=0)

hpt_chunk_file = []
for ci, beg_pi in enumerate(range(0, len(PS), chunk_size)):
    end_pi = min(beg_pi + chunk_size, len(PS))
    OPS = PS.iloc[beg_pi:end_pi]
    OPS.to_csv(f"{out_dir}/{ci}-pos.txt", index=False, sep="\t")

    hcf = open(f"{out_dir}/{ci}-hpth.txt", "wt")
    print("Family", "childId", "haplothreadId", "Haplothread", sep="\t", file=hcf)
    hpt_chunk_file.append((hcf, beg_pi, end_pi))


if in_hpth_file.endswith(".gz"):
    F = gzip.open(in_hpth_file, "rt")
else:
    F = open(in_hpth_file, "rt")

hcs = F.readline().strip("\n\r").split("\t")
assert len(hcs) == 4
for li, line in enumerate(F):
    if (li%100) == 0:
        print(f"Loading line {li}...")
    cs = line.strip("\n\r").split("\t")
    assert len(hcs) == len(cs)

    fm_id, p_id, hp_id, hp_s = cs
    assert len(hp_s) == len(PS)
    for hcf, beg_pi, end_pi in hpt_chunk_file:
        print(fm_id, p_id, hp_id, hp_s[beg_pi:end_pi], sep="\t", file=hcf)

F.close()
for hcf, _, _ in hpt_chunk_file:
    hcf.close()
