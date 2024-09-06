import itertools
import pandas as pd
import subprocess

import itertools
import pandas as pd
import subprocess
import random


ISOLATES = [i for i in open("/home/mahmadi/mpox_seqs/isolates_B1.txt").read().split('\n') if len(i) >0]

# Generate 4 random proportion combinations
isolate_combinations = list(itertools.combinations(ISOLATES, 2))
primer_bed_file = "/home/mahmadi/mpox_seqs/primer.bed"
for isolate_combination in isolate_combinations:
    proportion1 = round(random.uniform(0, 1), 2)
    proportion2 = round(1.0 - proportion1, 2)
    file1_path = "/home/mahmadi/mpox_seqs/sequences/"+ isolate_combination[0] + ".fasta"
    file2_path = "/home/mahmadi/mpox_seqs/sequences/" + isolate_combination[1] + ".fasta"
    output_path = "/home/mahmadi/mpox_seqs/simulations/combinations/B1/" + isolate_combination[0] + "_" + str(proportion1) + "_" +  isolate_combination[1] + "_" + str(proportion2)
    command = [
            "mixamp","simulate-proportions",
            f"{file1_path},{file2_path}",
            primer_bed_file,
            "--proportions", f"{proportion1},{proportion2}",
            "--outdir", output_path
        ]
        # Execute the command
    subprocess.run(command, check=True)


