#!/usr/bin/env python

from haplo_tools import load_haplothread_file
import pandas as pd
import sys

in_data_dir = "/mnt/wigclust8/data/safe/autism/HT/AGRE_WG38_AGRE859"
coll = "AGRE859"
chrom = "chr18"
chunk_size = 100000
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
RD = load_haplothread_file(in_hpth_file)

for ci, beg_pi in enumerate(range(0, len(PS), chunk_size)):
    end_pi = min(beg_pi + chunk_size, len(PS))
    OPS = PS.iloc[beg_pi:end_pi]
    OPS.to_csv(f"{out_dir}/{ci}-pos.txt", index=False, sep="\t")

    # with gzip.open(f"{out_dir}/{coll}-{chrom}-{ci}-hpth.txt.gz", "wt") as OF:
    with open(f"{out_dir}/{ci}-hpth.txt", "wt") as OF:
        print(f"Writing chunk {ci}...")
        print("Family", "childId", "haplothreadId", "Haplothread", sep="\t", file=OF)
        for fm_id, p_id, hp_id, hp_s in RD:
            assert len(hp_s) == len(PS)
            print(fm_id, p_id, hp_id, hp_s[beg_pi:end_pi], sep="\t", file=OF)