from functools import cache
import gzip

@cache
def load_haplothread_file(file_name: str) -> list[tuple[str,str,str,str]]:
    if file_name.endswith(".gz"):
        F = gzip.open(file_name, "rt")
    else:
        F = open(file_name, "rt")

    RD = []
    hcs = F.readline().strip("\n\r").split("\t")
    assert len(hcs) == 4
    for li, line in enumerate(F):
        if (li%100) == 0:
            print(f"Loading line {li}...")
        cs = line.strip("\n\r").split("\t")
        assert len(hcs) == len(cs)

        fmId, pId, hId, hps = cs
        RD.append((fmId, pId, hId, hps))

    F.close()
    return RD