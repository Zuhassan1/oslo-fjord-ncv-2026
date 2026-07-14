import re
from collections import defaultdict

SAMPLES = ["ERR14872774","ERR14872775","ERR14872776","ERR14872777","ERR14872778",
           "ERR14872779","ERR14872780","ERR14872781","ERR14872782"]
SAMPLE_META = {
    "ERR14872774": ("Operastranda", "T1", "flood"),
    "ERR14872775": ("Operastranda", "T2", "non-flood"),
    "ERR14872776": ("Operastranda", "T3", "non-flood"),
    "ERR14872777": ("Operastranda", "T4", "non-flood"),
    "ERR14872778": ("Operastranda", "T5", "non-flood"),
    "ERR14872779": ("Operastranda", "T6", "non-flood"),
    "ERR14872780": ("Operastranda", "T7", "non-flood"),
    "ERR14872781": ("Huk", "control", "control"),
    "ERR14872782": ("Huk", "control", "control"),
}

# Excluded: TIGTOG "Not_GV" with zero marker-gene support
EXCLUDE = {
    "ERR14872777_NODE_117_length_50263_cov_15.779637",
    "ERR14872777_NODE_216_length_35230_cov_18.204805",
}

with open("abundance_matrix_RPM.tsv") as f:
    header = f.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    totals_all = defaultdict(float)
    totals_filt = defaultdict(float)
    for line in f:
        fields = line.rstrip("\n").split("\t")
        cand = fields[0]
        vals = dict(zip(sample_cols, [float(x) for x in fields[1:]]))
        for col, v in vals.items():
            totals_all[col] += v
            if cand not in EXCLUDE:
                totals_filt[col] += v

print(f"{'Sample':<14}{'Site':<14}{'Category':<11}{'Total RPM (26)':>16}{'Total RPM (24, excl. NODE_117/216)':>36}{'%% change':>12}")
for s in SAMPLES:
    site, tp, cat = SAMPLE_META[s]
    a = totals_all[f"{s}_RPM"]
    b = totals_filt[f"{s}_RPM"]
    pct = ((b - a) / a * 100) if a > 0 else 0
    print(f"{s:<14}{site:<14}{cat:<11}{a:>16.2f}{b:>36.2f}{pct:>11.1f}%")
