from pathlib import Path
from snakeobjects import Project, ObjectGraph
import subprocess
import math

def run(proj: Project, OG: ObjectGraph):
    source_dir = Path(proj.parameters['source_dir'])
    collection = proj.parameters['collection']
    bin_size = proj.parameters['positions_per_bin']

    if "chromosomes" in proj.parameters:
        chromosomes = proj.parameters['chromosomes'].split(",")
    else:
        chromosomes = [f'chr{c}' for c in range(1,23)] # + ["chrX"]

    for chrom in chromosomes:
        in_chrom_o = OG.add('in_chromosome', chrom, {"chrom":chrom})
        pos_file = f"{source_dir}/{collection}-{chrom}-pos.txt"

        result = subprocess.run(['wc', '-l', pos_file], capture_output=True, text=True, check=True)
        line_count_str = result.stdout.strip().split(' ')[0]
        number_of_positions = int(line_count_str)-1
        n_bins = math.ceil(number_of_positions / bin_size)

        for chunk in range(n_bins):
            OG.add('chunk', f"{chrom}-{chunk}",{"chrom": chrom, "chunk": chunk}, [in_chrom_o])
