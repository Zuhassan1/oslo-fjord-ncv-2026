#!/usr/bin/env python3
import sys
from collections import defaultdict

NODES_DMP = "/hdd/zoh/databases/diamond_nr/nodes.dmp"
KRAKEN_OUT = "/hdd/zoh/oslo_virus_project/giant_virus_screening/candidate_catalog/kraken2_screening_output.tsv"
DIAMOND_OUT = "/hdd/zoh/oslo_virus_project/giant_virus_screening/candidate_catalog/diamond_screening_results.tsv"
PROTEIN_FAA = "/hdd/zoh/oslo_virus_project/giant_virus_screening/candidate_catalog/candidate_proteins.faa"

BACTERIA, ARCHAEA, EUKARYOTA, VIRUSES = "2", "2157", "2759", "10239"
CELLULAR_ROOTS = {BACTERIA, ARCHAEA, EUKARYOTA, "131567"}  # last = "cellular organisms", ambiguous-but-not-viral

print("Loading taxonomy nodes...", file=sys.stderr)
parent = {}
with open(NODES_DMP) as f:
    for line in f:
        parts = line.split("\t|\t")
        parent[parts[0].strip()] = parts[1].strip()
print(f"Loaded {len(parent)} nodes.", file=sys.stderr)

_cache = {}
def classify(taxid):
    if taxid in ("0", "A", ""):
        return "unclassified"
    if taxid in _cache:
        return _cache[taxid]
    seen, t, result = set(), taxid, "other"
    while t in parent and t not in seen:
        if t in CELLULAR_ROOTS:
            result = "cellular"; break
        if t == VIRUSES:
            result = "viral"; break
        seen.add(t)
        nxt = parent[t]
        if nxt == t:
            break
        t = nxt
    _cache[taxid] = result
    return result

print("Parsing Kraken2 output...", file=sys.stderr)
kraken_stats = {}
with open(KRAKEN_OUT) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        seqid, lca_string = fields[1], fields[4]
        counts = defaultdict(int)
        for token in lca_string.split():
            if ":" not in token:
                continue
            tid, cnt = token.rsplit(":", 1)
            try:
                cnt = int(cnt)
            except ValueError:
                continue
            counts[classify(tid)] += cnt
        kraken_stats[seqid] = dict(counts)
        kraken_stats[seqid]["total"] = sum(counts.values())

print("Counting proteins per contig...", file=sys.stderr)
protein_totals = defaultdict(int)
protein_to_contig = {}
with open(PROTEIN_FAA) as f:
    for line in f:
        if line.startswith(">"):
            pid = line[1:].split()[0]
            contig = pid.rsplit("_", 1)[0]
            protein_totals[contig] += 1
            protein_to_contig[pid] = contig
print(f"Total proteins counted: {sum(protein_totals.values())} (expect 1780)", file=sys.stderr)

print("Parsing DIAMOND output...", file=sys.stderr)
diamond_cellular = defaultdict(int)
diamond_hits = defaultdict(int)
total_diamond_lines = 0
with open(DIAMOND_OUT) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue
        total_diamond_lines += 1
        qseqid, staxids_field = fields[0], fields[6]
        first_taxid = staxids_field.split(";")[0].strip()
        contig = protein_to_contig.get(qseqid, qseqid.rsplit("_", 1)[0])
        diamond_hits[contig] += 1
        if classify(first_taxid) == "cellular":
            diamond_cellular[contig] += 1
print(f"Total DIAMOND hit lines: {total_diamond_lines} (expect 433)", file=sys.stderr)

all_contigs = sorted(set(list(protein_totals.keys()) + list(kraken_stats.keys())))
print(f"\nTotal contigs found: {len(all_contigs)} (expect 30)\n", file=sys.stderr)

print(f"{'contig':<58}{'kraken_%cell':>13}{'diamond_%cell':>15}{'n_prot':>8}{'n_hits':>8}")
results = []
for c in all_contigs:
    ks = kraken_stats.get(c)
    kp = (ks["cellular"] / ks["total"] * 100) if ks and ks.get("total", 0) > 0 else None
    n_prot = protein_totals.get(c, 0)
    n_hit = diamond_hits.get(c, 0)
    n_cell = diamond_cellular.get(c, 0)
    dp = (n_cell / n_prot * 100) if n_prot > 0 else None
    results.append((c, kp, dp, n_prot, n_hit))
    print(f"{c:<58}{(f'{kp:.2f}' if kp is not None else 'NA'):>13}{(f'{dp:.2f}' if dp is not None else 'NA'):>15}{n_prot:>8}{n_hit:>8}")

print("\n=== PASS/FAIL: reference paper thresholds (Kraken<=60% AND Diamond<=60%) ===")
n_pass = 0
for c, kp, dp, n_prot, n_hit in results:
    kp_pass = kp is not None and kp <= 60
    dp_pass = dp is not None and dp <= 60
    overall = kp_pass and dp_pass
    n_pass += overall
    print(f"{c:<58} Kraken<=60: {kp_pass!s:<6} Diamond<=60: {dp_pass!s:<6} PASS: {overall}")
print(f"\nCandidates passing BOTH filters: {n_pass} / {len(results)}")
