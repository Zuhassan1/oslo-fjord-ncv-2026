import re
from collections import defaultdict

SAMPLES = ["ERR14872774","ERR14872775","ERR14872776","ERR14872777","ERR14872778",
           "ERR14872779","ERR14872780","ERR14872781","ERR14872782"]
SAMPLE_META = {
    "ERR14872774": "Operastranda-T1-flood", "ERR14872775": "Operastranda-T2",
    "ERR14872776": "Operastranda-T3", "ERR14872777": "Operastranda-T4",
    "ERR14872778": "Operastranda-T5", "ERR14872779": "Operastranda-T6",
    "ERR14872780": "Operastranda-T7", "ERR14872781": "Huk-control-1",
    "ERR14872782": "Huk-control-2",
}

with open("abundance_matrix_RPM.tsv") as f:
    header = f.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    rows = []
    for line in f:
        fields = line.rstrip("\n").split("\t")
        cand = fields[0]
        origin = re.match(r'^(ERR\d+)_', cand).group(1)
        vals = dict(zip(sample_cols, [float(x) for x in fields[1:]]))
        rows.append((cand, origin, vals))

print("Candidate count by origin sample (confirms the imbalance):")
origin_counts = defaultdict(int)
for cand, origin, _ in rows:
    origin_counts[origin] += 1
for s in SAMPLES:
    print(f"  {s} ({SAMPLE_META[s]}): {origin_counts[s]} candidates contributed")

print(f"\n{'Candidate':<58}{'Origin':<14}", "\t".join(SAMPLE_META[s] for s in SAMPLES))
for cand, origin, vals in rows:
    marker = "*" if origin else ""
    line_vals = []
    for s in SAMPLES:
        v = vals[f"{s}_RPM"]
        tag = "[HOME]" if s == origin else ""
        line_vals.append(f"{v:.2f}{tag}")
    print(f"{cand:<58}{origin:<14}" + "\t".join(line_vals))
