#!/usr/bin/env python

import sys
import gzip

out_dir = "./individuals"
vcfs = [
    "/mnt/wigtop2/data/unsafe/iossifov/AGRE_to_vcf/chunk/chr9-0/chunk.vcf",
    "/mnt/wigtop2/data/unsafe/iossifov/AGRE_to_vcf/chunk/chr9-1/chunk.vcf",
    "/mnt/wigtop2/data/unsafe/iossifov/AGRE_to_vcf/chunk/chr9-2/chunk.vcf"
]

assert len(sys.argv) != 1
out_dir = sys.argv[1]
vcfs = sys.argv[2:]

pos_f = open(f"{out_dir}/pos.txt", "wt")
print("CHROM", "POS", "REF", "ALT", sep="\t", file=pos_f)
samples = []
sample_files = []
for vcf in vcfs:
    print(f"Processing {vcf}...")
    with gzip.open(vcf,"rt") as V:
        line = ""
        for li,line in enumerate(V):
            if line.startswith("#CHROM"):
                break
            assert li < 10
        assert line.startswith("#CHROM")
        hcs = line.strip("\n\r").split("\t")
        v_samples = hcs[9:]
        assert not samples or samples == v_samples
        if not samples:
            samples = v_samples
            sample_files = [gzip.open(f"{out_dir}/{sample}-gens.txt.gz", "wt") for sample in samples]
        for li, line in enumerate(V, li+1):
            if (li % 1000) == 0:
                print(f"   Working on line {li}")
            cs = line.strip("\n\r").split("\t")
            print(cs[0], cs[1], cs[3], cs[4], sep="\t", file=pos_f)
            assert len(cs) == len(hcs)
            gens = cs[9:]
            for sample_file, gen in zip(sample_files, gens):
                print(gen, file=sample_file)

for sample_file in sample_files:
    sample_file.close()
pos_f.close()
