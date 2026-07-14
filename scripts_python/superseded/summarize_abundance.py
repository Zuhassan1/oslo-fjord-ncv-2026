import re
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

with open("abundance_matrix_RPM.tsv") as f:
    header = f.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]  # e.g. ERR14872774_RPM, ...

    totals = defaultdict(float)
    per_candidate = {}
    for line in f:
        fields = line.rstrip("\n").split("\t")
        cand = fields[0]
        vals = [float(x) for x in fields[1:]]
        per_candidate[cand] = dict(zip(sample_cols, vals))
        for col, v in zip(sample_cols, vals):
            totals[col] += v

print(f"{'Sample':<14}{'Site':<14}{'Timepoint':<11}{'Date':<12}{'Category':<11}{'Total RPM (all 26 candidates)':>30}")
for s in SAMPLES:
    site, tp, date, cat = SAMPLE_META[s]
    total_rpm = totals[f"{s}_RPM"]
    print(f"{s:<14}{site:<14}{tp:<11}{date:<12}{cat:<11}{total_rpm:>30.2f}")

print("\n=== Top 3 contributing candidates per sample (by RPM) ===")
for s in SAMPLES:
    col = f"{s}_RPM"
    ranked = sorted(per_candidate.items(), key=lambda kv: kv[1][col], reverse=True)[:3]
    print(f"\n{s}:")
    for cand, vals in ranked:
        if vals[col] > 0:
            print(f"  {cand}: {vals[col]:.3f} RPM")
