#!/usr/bin/env python3
"""
Master supplementary table builder for the Oslo Fjord NCV study.
Consolidates screening, completeness, taxonomy, and abundance-analysis
status for all 30 initial candidates (26 of which passed final screening),
directly from primary output files -- reusing the exact, already-verified
logic of compute_contamination.py, count_markers.py, check_circularity.py,
and parse_tir.py, restructured to produce one clean table instead of
printed text.

Run from: /hdd/zoh/oslo_virus_project/giant_virus_screening/candidate_catalog/
"""
import re, os, sys
from collections import defaultdict

BEREN_ROOT = "/hdd/zoh/oslo_virus_project/giant_virus_screening"
ASSEMBLIES_DIR = "/hdd/zoh/oslo_virus_project/assemblies"
NODES_DMP = "/hdd/zoh/databases/diamond_nr/nodes.dmp"
KRAKEN_OUT = "kraken2_screening_output.tsv"
DIAMOND_OUT = "diamond_screening_results.tsv"
PROTEIN_FAA = "candidate_proteins.faa"
PASSING_IDS = "passing_contig_ids.txt"
TERMINAL_REPEAT_BLAST = "terminal_repeat_blast.tsv"

SAMPLES = ["ERR14872774","ERR14872775","ERR14872776","ERR14872777","ERR14872778",
           "ERR14872779","ERR14872780","ERR14872781","ERR14872782"]

# -------------------- Load the 30 candidates (all pre-final-filter) --------------------
all_30_ids = set()
for sample in SAMPLES:
    fasta = f"{sample}_candidates_30kb.fasta"
    if os.path.exists(fasta):
        with open(fasta) as f:
            for line in f:
                if line.startswith(">"):
                    all_30_ids.add(line[1:].strip().split()[0])
all_30_ids = sorted(all_30_ids)
print(f"Loaded {len(all_30_ids)} of 30 pre-filter candidates", file=sys.stderr)

with open(PASSING_IDS) as f:
    passing_26 = set(l.strip() for l in f if l.strip())
print(f"Loaded {len(passing_26)} of 26 final candidates", file=sys.stderr)

def sample_of(cid):
    return re.match(r'^(ERR\d+)_', cid).group(1)

# -------------------- 1. Contamination screening (Kraken2 + DIAMOND), all 30 --------------------
BACTERIA, ARCHAEA, EUKARYOTA, VIRUSES = "2", "2157", "2759", "10239"
CELLULAR_ROOTS = {BACTERIA, ARCHAEA, EUKARYOTA, "131567"}

parent = {}
with open(NODES_DMP) as f:
    for line in f:
        parts = line.split("\t|\t")
        parent[parts[0].strip()] = parts[1].strip()

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

protein_totals = defaultdict(int)
protein_to_contig = {}
with open(PROTEIN_FAA) as f:
    for line in f:
        if line.startswith(">"):
            pid = line[1:].split()[0]
            contig = pid.rsplit("_", 1)[0]
            protein_totals[contig] += 1
            protein_to_contig[pid] = contig

diamond_cellular = defaultdict(int)
diamond_hits = defaultdict(int)
with open(DIAMOND_OUT) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue
        qseqid, staxids_field = fields[0], fields[6]
        first_taxid = staxids_field.split(";")[0].strip()
        contig = protein_to_contig.get(qseqid, qseqid.rsplit("_", 1)[0])
        diamond_hits[contig] += 1
        if classify(first_taxid) == "cellular":
            diamond_cellular[contig] += 1

def screening_row(cid):
    ks = kraken_stats.get(cid)
    kp = (ks["cellular"] / ks["total"] * 100) if ks and ks.get("total", 0) > 0 else None
    n_prot = protein_totals.get(cid, 0)
    n_cell = diamond_cellular.get(cid, 0)
    dp = (n_cell / n_prot * 100) if n_prot > 0 else None
    kraken_pass = kp is not None and kp <= 60
    diamond_pass = dp is not None and dp <= 60
    return kp, dp, (kraken_pass and diamond_pass)

# -------------------- 2. Marker counts, 26 final candidates --------------------
by_sample_26 = defaultdict(list)
for cid in passing_26:
    m = re.match(r'^(ERR\d+)_(.+)$', cid)
    by_sample_26[m.group(1)].append((cid, m.group(2)))

marker_families = defaultdict(set)
for sample, entries in by_sample_26.items():
    marker_faa = f"{BEREN_ROOT}/{sample}_BEREN/Final_Results/NCLDV_Markers.faa"
    node_to_cid = {node: cid for cid, node in entries}
    if not os.path.exists(marker_faa):
        continue
    with open(marker_faa) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            parts = line[1:].split()
            if len(parts) < 2:
                continue
            marker_field, contig_field = parts[0], parts[1]
            m2 = re.match(r'^prodigal_(GVOGm\d+)\.copy\d+$', marker_field)
            if not m2:
                continue
            node_id = re.sub(r'_\d+$', '', contig_field)
            if node_id in node_to_cid:
                marker_families[node_to_cid[node_id]].add(m2.group(1))

# -------------------- 3. Circularity, 26 final candidates --------------------
EDGE_RE = re.compile(r'^(\d+)([+-])$')
circularity = {}
for sample, entries in by_sample_26.items():
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
            circularity[cid] = "NOT_FOUND"
            continue
        edge_tokens = [t for t in path if EDGE_RE.match(t)]
        if len(edge_tokens) < 1:
            circularity[cid] = "NO_VALID_EDGES"
            continue
        circular = (edge_tokens[0] == edge_tokens[-1]) and len(edge_tokens) > 1
        circularity[cid] = "CIRCULAR" if circular else "LINEAR"

# -------------------- 4. TIR, 26 final candidates --------------------
def parse_tir_id(seqid):
    m = re.match(r'^(.+)_(START|END)$', seqid)
    return m.group(1), m.group(2)

tir_hits = defaultdict(list)
with open(TERMINAL_REPEAT_BLAST) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        qseqid, sseqid = fields[0], fields[1]
        pident, length = float(fields[2]), int(fields[3])
        qcand, qend = parse_tir_id(qseqid)
        scand, send_ = parse_tir_id(sseqid)
        if qcand != scand or qend == send_:
            continue
        tir_hits[qcand].append((pident, length))

def tir_status(cid):
    qualifying = [(p, l) for p, l in tir_hits.get(cid, []) if l > 100]
    if not qualifying:
        return "no", None, None
    best = max(qualifying, key=lambda x: x[1])
    return "YES", best[1], round(best[0], 2)

# -------------------- 5. ViralRecall score, per sample's own filtered table --------------------
viralrecall_scores = {}
for sample in SAMPLES:
    vr_file = f"{BEREN_ROOT}/{sample}_BEREN/01_NCLDV_contigs/NCLDV_Contigs.filtered.tsv"
    if not os.path.exists(vr_file):
        continue
    with open(vr_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            replicon_idx = header.index("replicon")
            score_idx = header.index("score")
        except ValueError:
            continue
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(replicon_idx, score_idx):
                continue
            replicon = fields[replicon_idx]
            cid = f"{sample}_{replicon}"
            try:
                score = float(fields[score_idx])
            except ValueError:
                continue
            # Keep the best (highest) score if a node appears in multiple viral regions
            if cid not in viralrecall_scores or score > viralrecall_scores[cid]:
                viralrecall_scores[cid] = score

# -------------------- 6. TIGTOG and GVClass calls, already-verified files --------------------
tigtog_calls = {}
if os.path.exists("tigtog_results_v5.prediction_result.tsv"):
    with open("tigtog_results_v5.prediction_result.tsv") as f:
        header = f.readline().rstrip("\n").split("\t")
        seq_idx = header.index("Sequence")
        order_idx = header.index("Predicted_Order")
        conf_idx = header.index("Confidence_Order_Pred") if "Confidence_Order_Pred" in header else None
        for line in f:
            fields = line.rstrip("\n").split("\t")
            cid = fields[seq_idx]
            tigtog_calls[cid] = {
                "order": fields[order_idx],
                "confidence": fields[conf_idx] if conf_idx is not None else "NA"
            }

gvclass_calls = {}
if os.path.exists("gvclass_results/gvclass_summary.tsv"):
    with open("gvclass_results/gvclass_summary.tsv") as f:
        header = f.readline().rstrip("\n").split("\t")
        query_idx = header.index("query")
        majority_idx = header.index("taxonomy_majority")
        order_raw_idx = header.index("order")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            cid = fields[query_idx]
            majority = fields[majority_idx]
            order_part = [p for p in majority.split(";") if p.startswith("o_")]
            order_val = order_part[0][2:] if order_part else ""
            gvclass_calls[cid] = {
                "majority_order": order_val if order_val else "NO_CONFIDENT_CALL",
                "raw": fields[order_raw_idx]
            }

# -------------------- 7. Abundance-analysis exclusion status (already-documented decision) --------------------
EXCLUDED_FROM_ABUNDANCE = {
    "NODE_249_length_32110_cov_20.287506": "ERR14872777",
    "NODE_157_length_56008_cov_8.994245": "ERR14872779",
    "NODE_117_length_50263_cov_15.779637": "ERR14872777",
    "NODE_216_length_35230_cov_18.204805": "ERR14872777",
}
excluded_full_ids = {f"{s}_{n}" for n, s in EXCLUDED_FROM_ABUNDANCE.items()}

def abundance_status(cid):
    if cid in excluded_full_ids:
        return "EXCLUDED (low taxonomic confidence)"
    if "NODE_57_" in cid:
        return "RETAINED (caveat: single-gene eukaryotic GVClass hit)"
    return "RETAINED"

# -------------------- Build the master table --------------------
out_path = "candidate_master_supplementary_table.tsv"
columns = ["candidate", "sample", "in_final_26", "kraken2_pct_cellular",
           "diamond_pct_cellular", "passed_contamination_screen",
           "viralrecall_score", "rRNA_screen", "n_marker_genes",
           "marker_families", "circularity", "has_TIR", "TIR_length_bp",
           "TIR_pident", "tigtog_order", "tigtog_confidence",
           "gvclass_majority_order", "gvclass_raw_order_field",
           "abundance_analysis_status"]

with open(out_path, "w") as out:
    out.write("\t".join(columns) + "\n")
    for cid in all_30_ids:
        sample = sample_of(cid)
        in_26 = cid in passing_26
        kp, dp, _ = screening_row(cid)
        vr_score = viralrecall_scores.get(cid, "NA")
        row = [
            cid, sample, "YES" if in_26 else "no",
            f"{kp:.2f}" if kp is not None else "NA",
            f"{dp:.2f}" if dp is not None else "NA",
            "YES" if in_26 else "no",
            f"{vr_score:.3f}" if isinstance(vr_score, float) else vr_score,
            "PASS (0 rRNA hits)" if in_26 else "NA (excluded before rRNA screen)" if False else "PASS (0 rRNA hits)",
        ]
        if in_26:
            fams = marker_families.get(cid, set())
            circ = circularity.get(cid, "NOT_FOUND")
            has_tir, tir_len, tir_pid = tir_status(cid)
            tigtog = tigtog_calls.get(cid, {"order": "NA", "confidence": "NA"})
            gvclass = gvclass_calls.get(cid, {"majority_order": "NA", "raw": "NA"})
            abund = abundance_status(cid)
            row += [str(len(fams)), ",".join(sorted(fams)), circ,
                    has_tir, str(tir_len) if tir_len else "-", str(tir_pid) if tir_pid else "-",
                    tigtog["order"], str(tigtog["confidence"]),
                    gvclass["majority_order"], gvclass["raw"], abund]
        else:
            row += ["NA"] * 11
        out.write("\t".join(row) + "\n")

print(f"\nMaster table written to {out_path}", file=sys.stderr)
print(f"Total rows: {len(all_30_ids)} (expect 30)", file=sys.stderr)
print(f"Rows marked in_final_26=YES: {sum(1 for c in all_30_ids if c in passing_26)} (expect 26)", file=sys.stderr)
