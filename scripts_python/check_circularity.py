#!/usr/bin/env python3
import re, sys
from collections import defaultdict

ASSEMBLIES_DIR = "/hdd/zoh/oslo_virus_project/assemblies"
PASSING_IDS = "/hdd/zoh/oslo_virus_project/giant_virus_screening/candidate_catalog/passing_contig_ids.txt"
EDGE_RE = re.compile(r'^(\d+)([+-])$')

with open(PASSING_IDS) as f:
    candidate_ids = [line.strip() for line in f if line.strip()]
print(f"Loaded {len(candidate_ids)} candidate IDs (expect 26)", file=sys.stderr)

by_sample = defaultdict(list)
for cid in candidate_ids:
    m = re.match(r'^(ERR\d+)_(.+)$', cid)
    if not m:
        print(f"WARNING: could not parse sample from {cid}", file=sys.stderr)
        continue
    by_sample[m.group(1)].append((cid, m.group(2)))

results = {}

for sample, entries in by_sample.items():
    paths_file = f"{ASSEMBLIES_DIR}/{sample}_assembly/scaffolds.paths"
    wanted = {node for _, node in entries}
    found = {}
    header, is_revcomp, tokens = None, False, []

    def flush():
        if header is not None and not is_revcomp and header in wanted:
            found[header] = list(tokens)

    with open(paths_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("NODE_"):
                flush()
                is_revcomp = line.endswith("'")
                header = line[:-1] if is_revcomp else line
                tokens = []
            else:
                tokens.extend(t for t in line.split(",") if t)
        flush()

    for cid, node in entries:
        path = found.get(node)
        if path is None:
            results[cid] = ("NOT_FOUND", 0, None, None, 0)
            continue
        edge_tokens = [t for t in path if EDGE_RE.match(t)]
        skipped = len(path) - len(edge_tokens)
        if len(edge_tokens) < 1:
            results[cid] = ("NO_VALID_EDGES", len(path), None, None, skipped)
            continue
        first, last = edge_tokens[0], edge_tokens[-1]
        circular = (first == last) and len(edge_tokens) > 1
        results[cid] = ("CIRCULAR" if circular else "LINEAR", len(edge_tokens), first, last, skipped)

print(f"\n{'contig':<58}{'status':<14}{'n_edges':>8}{'first':>10}{'last':>10}{'skipped':>9}")
for cid in candidate_ids:
    status, n_edges, first, last, skipped = results.get(cid, ("MISSING", 0, None, None, 0))
    print(f"{cid:<58}{status:<14}{n_edges:>8}{str(first):>10}{str(last):>10}{skipped:>9}")

n_circ = sum(1 for v in results.values() if v[0] == "CIRCULAR")
n_lin = sum(1 for v in results.values() if v[0] == "LINEAR")
n_other = len(results) - n_circ - n_lin
print(f"\nTotal processed: {len(results)} (expect 26)")
print(f"Circular: {n_circ}  |  Linear: {n_lin}  |  Other/unresolved: {n_other}")
