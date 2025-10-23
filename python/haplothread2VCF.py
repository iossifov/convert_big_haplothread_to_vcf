#!/usr/bin/env python

import sys
import os
from collections import defaultdict
import pandas as pd
import gzip

from haplo_tools import load_haplothread_file


in_data_dir = "/mnt/wigclust8/data/safe/autism/HT/AGRE_WG38_AGRE859"
coll = "AGRE859"
chrom = "chr18"

meta_data_file = f"{in_data_dir}/metaData.txt"
in_pos_file = f"{in_data_dir}/{coll}-{chrom}-pos.txt"
in_hpth_file = f"{in_data_dir}/{coll}-{chrom}-hpth.txt.gz"
out_file = f"out-vcfs/{coll}_{chrom}.vcf"

assert len(sys.argv) in [1, 6]
if len(sys.argv) == 6:
    meta_data_file, in_pos_file, in_hpth_file, chrom, out_file = sys.argv[1:]

print(f"Working on {meta_data_file=}, {in_pos_file=}, {in_hpth_file=}, {chrom=} {out_file=}...")

PS = pd.read_table(in_pos_file, sep='\t', header=0)
RD = load_haplothread_file(in_hpth_file)

hap_idx = {(p_id, h_id):idx for idx, (_, p_id, h_id, _) in enumerate(RD)}

MT = pd.read_table(meta_data_file, sep='\t', header=0)
person_tuple = {t.indId:t for t in MT.itertuples(index=False)}

child_haps_idx: dict[tuple[str, str, str], dict[str, int]] = defaultdict(dict)
for h_idx, (fm_id, p_id, h_id, _) in enumerate(RD):
    m_id, f_id = fm_id.split("_")
    child_haps_idx[m_id, f_id, p_id][h_id] = h_idx
for hd in child_haps_idx.values():
    assert len(hd) == 4
    assert set(hd.keys()) == set(['MT', 'MNT', 'FT', 'FNT'])

person_hap_idxs: dict[str, list[tuple[int,int]]] = defaultdict(list)
for (m_id, f_id, p_id), hd in child_haps_idx.items():
    person_hap_idxs[p_id].append((hd['MT'], hd['FT'], f"child"))
    person_hap_idxs[m_id].append((hd['MT'], hd['MNT'], f"mother of {p_id}"))
    person_hap_idxs[f_id].append((hd['FT'], hd['FNT'], f"father of {p_id}"))


    # atts = eval(person_tuple[p_id].atts)
    # a_mom_id = atts['mom']
    # a_dad_id = atts['dad']
    # if (person_tuple[p_id].dad != f_id) or \
    #    (person_tuple[p_id].mom != m_id):
    #     print("TTTTT", p_id, h_id, f_id, m_id, person_tuple[p_id].dad, person_tuple[p_id].mom, a_dad_id, a_mom_id)
    # if (a_dad_id != f_id) or \
    #    (a_mom_id != m_id):
    #     print("GGGGG", p_id, h_id, f_id, m_id, person_tuple[p_id].dad, person_tuple[p_id].mom, a_dad_id, a_mom_id)

chrLens = {
    "chr1":   248_956_422,
    "chr2":   242_193_529,
    "chr3":   198_295_559,
    "chr4":   190_214_555,
    "chr5":   181_538_259,
    "chr6":   170_805_979,
    "chr7":   159_345_973,
    "chr8":   145_138_636,
    "chr9":   138_394_717,
    "chr10":  133_797_422,
    "chr11":  135_086_622,
    "chr12":  133_275_309,
    "chr13":  114_364_328,
    "chr14":  107_043_718,
    "chr15":  101_991_189,
    "chr16":   90_338_345,
    "chr17":   83_257_441,
    "chr18":   80_373_285,
    "chr19":   58_617_616,
    "chr20":   64_444_167,
    "chr21":   46_709_983,
    "chr22":   50_818_468,
    "chrX":   156_040_895
}



ambiguous = {
    'M': set('AC'),
    'R': set('AG'),
    'W': set('AT'),
    'S': set('CG'),
    'Y': set('CT'),
    'K': set('GT')
}

def get_alleles(a):
    if a in 'codu':
        return set([])
    if a in 'ACGT':
        return set([a])
    if a in ambiguous:
        return ambiguous[a]
    assert False


def toVCFGen(a, b, refA):
    if a in "ACGT":
        assert b in "ACGT"
        if a == b:
            if a == refA:
                return "0/0"
            else:
                return "1/1"
        else:
            return "0/1"
    else:
        assert a == b
        if a in 'codu':
            return "./."
        else:
            return "0/1"


def parGen(gens):
    if not gens:
        return "./."
    if len(gens) > 1:
        gens.remove('./.')
    gn, = gens
    return gn



samples = list(person_hap_idxs.keys())

os.makedirs(os.path.dirname(out_file), exist_ok=True)

# open gzipped output when file ends with .gz
if out_file.endswith('.gz'):
    OF = gzip.open(out_file, 'wt')
else:
    OF = open(out_file, 'w')


print('##fileformat=VCFv4.2', file=OF)
print('##fileDate=20251021', file=OF)
print('##source=haplothread2VCF.py', file=OF)
print('##reference=GRCh38-hg38', file=OF)
print(f'##contig=<ID={chrom},length={chrLens[chrom]}>', file=OF)
print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=OF)

print(*('#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT'.split(' ')), *samples,
      sep="\t", file=OF)





for pi, t in enumerate(PS.itertuples(index=False)):
    ref_allele = t.RefAllele

    pos_alleles = {a for c, a in zip([t.A, t.C, t.G, t.T], "ACGT") if c > 0}
    if not pos_alleles:
        print(f"WARNING: no alleles observed at {pi=}, {t=}")
        continue
    pos_alt_alleles = set(pos_alleles)
    if ref_allele in pos_alt_alleles:
        pos_alt_alleles.remove(ref_allele)
    if len(pos_alt_alleles) > 1:
        print(f"WARNING: multiple alternative alleles found at {pi=}, {t=}, {pos_alt_alleles=}")


    data_alleles = set()
    person_allele_set = {}
    for p_id in samples:
        hap_idx = person_hap_idxs[p_id]
        person_allele_sets = set([])
        for h1_idx, h2_idx, _ in hap_idx:
            alleles = get_alleles(RD[h1_idx][3][pi]) | get_alleles(RD[h2_idx][3][pi])
            assert len(alleles) < 3
            person_allele_sets.add(tuple(alleles))
        if len(person_allele_sets) == 2:
            if () in person_allele_sets:
                person_allele_sets.remove(())
            else:
                print(f"WARNING: weird person_allele_sets, {p_id=}, {pi=}, {person_allele_sets=}")
        assert len(person_allele_sets) > 0
        if len(person_allele_sets) > 1:
            person_alleles = sorted(person_allele_sets, key=lambda x: len(x))[-1]
            print(f"WARNING: more than one person alleles, {p_id=}, {pi=}, {person_allele_sets=}. Will use {person_alleles=}")
        else:
            person_alleles, = person_allele_sets
        data_alleles |= set(person_alleles)
        person_allele_set[p_id] = set(person_alleles)
    assert data_alleles.issubset(pos_alleles)
    if data_alleles != pos_alleles:
        print("DPDPDP", data_alleles, pos_alleles, ref_allele)
    data_alt_alleles = sorted(data_alleles - set([ref_allele]))
    data_alleles_d = {a: an for an, a in enumerate(data_alt_alleles,1)}
    data_alleles_d[ref_allele] = 0
    person_gens = {}
    for p_id, person_alleles in person_allele_set.items():
        allele_is = sorted([data_alleles_d[a] for a in person_alleles])
        if len(allele_is) == 0:
            person_gens[p_id] = "./."
            continue
        if len(allele_is) == 1:
            allele_is.append(allele_is[0])
        assert len(allele_is) == 2
        person_gens[p_id] = f"{allele_is[0]}/{allele_is[1]}"

    # if not data_alt_alleles:
    #     continue
    cs = [chrom, t.Position, ".", ref_allele, ",".join(data_alt_alleles),
        ".", ".", ".", "GT"] + [person_gens[p_id] for p_id in samples]
    print(*cs, sep="\t", file=OF)
    # assert len(data_alleles_d)<3

OF.close()
