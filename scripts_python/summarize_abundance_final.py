import re
from collections import defaultdict

SAMPLES = ["ERR14872774","ERR14872775","ERR14872776","ERR14872777","ERR14872778",
           "ERR14872779","ERR14872780","ERR14872781","ERR14872782"]
SAMPLE_META = {
    "ERR14872774": ("Operastranda", "flood"), "ERR14872775": ("Operastranda", "non-flood"),
    "ERR14872776": ("Operastranda", "non-flood"), "ERR14872777": ("Operastranda", "non-flood"),
    "ERR14872778": ("Operastranda", "non-flood"), "ERR14872779": ("Operastranda", "non-flood"),
    "ERR14872780": ("Operastranda", "non-flood"), "ERR14872781": ("Huk", "control"),
    "ERR14872782": ("Huk", "control"),
}

# Consistent with the taxonomy section's own documented conclusions:
# candidates with zero/near-zero marker support AND a Not_GV classifier call
EXCLUDE = {
    "ERR14872777_NODE_249_length_32110_cov_20.287506",   # strongest joint negative: blank GVClass + TIGTOG 1.00
    "ERR14872779_NODE_157_length_56008_cov_8.994245",    # strong joint negative: blank GVClass + TIGTOG 0.80
    "ERR14872777_NODE_117_length_50263_cov_15.779637",   # weak GVClass (0 markers) + TIGTOG 0.79; near-uniform cross-sample recruitment pattern
    "ERR14872777_NODE_216_length_35230_cov_18.204805",   # weak GVClass (0 markers) + TIGTOG 0.89
}
# Retained per earlier documented reasoning: NODE_333, NODE_155 (real multi-marker support
# despite TIGTOG disagreement), NODE_57 (passed primary filters, flagged with caveat only)

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

print(f"Excluded from this corrected view: {len(EXCLUDE)} candidates (see reasons above)")
print(f"Retained: {len(rows) - len(EXCLUDE)} of {len(rows)} candidates\n")

# Corrected totals
totals = defaultdict(float)
for cand, origin, vals in rows:
    if cand in EXCLUDE:
        continue
    for col, v in vals.items():
        totals[col] += v

print(f"{'Sample':<14}{'Site':<14}{'Category':<11}{'Corrected Total RPM (22 candidates)':>36}")
for s in SAMPLES:
    site, cat = SAMPLE_META[s]
    print(f"{s:<14}{site:<14}{cat:<11}{totals[f'{s}_RPM']:>36.2f}")

# Home-dominance ratio per retained candidate: HOME value / (best non-home value)
print("\n=== Home-dominance ratio per retained candidate (HOME RPM / best cross-sample RPM) ===")
print("(High ratio = tightly localized to one sample; low ratio = broadly distributed)")
for cand, origin, vals in rows:
    if cand in EXCLUDE:
        continue
    home_val = vals[f"{origin}_RPM"]
    others = [v for col, v in vals.items() if col != f"{origin}_RPM"]
    best_other = max(others) if others else 0
    ratio = (home_val / best_other) if best_other > 0 else float('inf')
    ratio_str = f"{ratio:.1f}x" if ratio != float('inf') else "inf (only in home sample)"
    print(f"  {cand:<58} home={home_val:.1f}  best_other={best_other:.1f}  ratio={ratio_str}")
