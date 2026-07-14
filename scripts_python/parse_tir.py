import re
from collections import defaultdict

def parse_id(seqid):
    m = re.match(r'^(.+)_(START|END)$', seqid)
    return m.group(1), m.group(2)

hits_by_candidate = defaultdict(list)
with open("terminal_repeat_blast.tsv") as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        qseqid, sseqid = fields[0], fields[1]
        pident, length = float(fields[2]), int(fields[3])
        qcand, qend = parse_id(qseqid)
        scand, send_ = parse_id(sseqid)
        if qcand != scand or qend == send_:
            continue
        hits_by_candidate[qcand].append((pident, length))

with open("passing_contig_ids.txt") as f:
    all_candidates = [l.strip() for l in f if l.strip()]

print(f"{'contig':<58}{'has_TIR':<10}{'best_len':>10}{'best_pident':>12}")
n_tir = 0
for c in all_candidates:
    qualifying = [(p, l) for p, l in hits_by_candidate.get(c, []) if l > 100]
    has_tir = len(qualifying) > 0
    n_tir += has_tir
    if qualifying:
        best = max(qualifying, key=lambda x: x[1])
        print(f"{c:<58}{'YES':<10}{best[1]:>10}{best[0]:>12.2f}")
    else:
        print(f"{c:<58}{'no':<10}{'-':>10}{'-':>12}")

print(f"\nCandidates with terminal inverted repeats (>100bp, minus-strand): {n_tir} / {len(all_candidates)}")
