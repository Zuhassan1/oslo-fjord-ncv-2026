import re, os
from collections import defaultdict

SAMPLES = ["ERR14872774","ERR14872775","ERR14872776","ERR14872777","ERR14872778",
           "ERR14872779","ERR14872780","ERR14872781","ERR14872782"]

SAMPLE_META = {
    "ERR14872774": ("Operastranda", "T1", "28-05-2024", "flood"),
    "ERR14872775": ("Operastranda", "T2", "04-06-2024", "non-flood"),
    "ERR14872776": ("Operastranda", "T3", "18-06-2024", "non-flood"),
    "ERR14872777": ("Operastranda", "T4", "02-07-2024", "non-flood"),
    "ERR14872778": ("Operastranda", "T5", "09-07-2024", "non-flood"),
    "ERR14872779": ("Operastranda", "T6", "06-08-2024", "non-flood"),
    "ERR14872780": ("Operastranda", "T7", "20-08-2024", "non-flood"),
    "ERR14872781": ("Huk", "control", "04-06-2024", "control"),
    "ERR14872782": ("Huk", "control", "02-07-2024", "control"),
}

ASSEMBLIES_DIR = "/hdd/zoh/oslo_virus_project/assemblies"

# Pull real total read counts from each sample's own assembly log
total_reads = {}
for s in SAMPLES:
    log_path = f"{ASSEMBLIES_DIR}/{s}_assembly.log"
    with open(log_path) as f:
        content = f.read()
    matches = re.findall(r'Total (\d+) reads processed', content)
    if not matches:
        raise ValueError(f"Could not find total read count for {s}")
    total_reads[s] = int(matches[0])
    print(f"{s}: {total_reads[s]:,} total reads (from assembly log)")

# Parse each sample's coverage.tsv
candidate_data = defaultdict(dict)  # candidate_data[candidate][sample] = (numreads, meandepth)
for s in SAMPLES:
    with open(f"{s}_coverage.tsv") as f:
        header = f.readline()
        for line in f:
            fields = line.rstrip("\n").split("\t")
            rname, numreads, meandepth = fields[0], int(fields[3]), float(fields[6])
            candidate_data[rname][s] = (numreads, meandepth)

candidates = sorted(candidate_data.keys())
print(f"\nTotal candidates found across coverage files: {len(candidates)} (expect 26)")

# Build normalized matrix: reads-per-million-total-reads (RPM), our own disclosed addition
with open("abundance_matrix_RPM.tsv", "w") as out:
    out.write("candidate\t" + "\t".join(f"{s}_RPM" for s in SAMPLES) + "\n")
    for cand in candidates:
        row = [cand]
        for s in SAMPLES:
            numreads, _ = candidate_data[cand].get(s, (0, 0.0))
            rpm = (numreads / total_reads[s]) * 1_000_000
            row.append(f"{rpm:.3f}")
        out.write("\t".join(row) + "\n")

with open("abundance_matrix_meandepth.tsv", "w") as out:
    out.write("candidate\t" + "\t".join(f"{s}_depth" for s in SAMPLES) + "\n")
    for cand in candidates:
        row = [cand]
        for s in SAMPLES:
            _, depth = candidate_data[cand].get(s, (0, 0.0))
            row.append(f"{depth:.4f}")
        out.write("\t".join(row) + "\n")

print("\nWrote abundance_matrix_RPM.tsv and abundance_matrix_meandepth.tsv")
