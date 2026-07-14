import re, os
from collections import defaultdict

BEREN_ROOT = "/hdd/zoh/oslo_virus_project/giant_virus_screening"
PASSING_IDS = "/hdd/zoh/oslo_virus_project/giant_virus_screening/candidate_catalog/passing_contig_ids.txt"

with open(PASSING_IDS) as f:
    candidate_ids = [l.strip() for l in f if l.strip()]
print(f"Loaded {len(candidate_ids)} candidate IDs (expect 26)")

by_sample = defaultdict(list)
for cid in candidate_ids:
    m = re.match(r'^(ERR\d+)_(.+)$', cid)
    sample, node = m.group(1), m.group(2)
    by_sample[sample].append((cid, node))

marker_families = defaultdict(set)

for sample, entries in by_sample.items():
    marker_faa = f"{BEREN_ROOT}/{sample}_BEREN/Final_Results/NCLDV_Markers.faa"
    node_to_cid = {node: cid for cid, node in entries}
    if not os.path.exists(marker_faa):
        print(f"WARNING: missing {marker_faa}")
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
            family = m2.group(1)
            node_id = re.sub(r'_\d+$', '', contig_field)
            if node_id in node_to_cid:
                marker_families[node_to_cid[node_id]].add(family)

print(f"\n{'contig':<58}{'n_markers':>11}   families")
for cid in candidate_ids:
    fams = marker_families.get(cid, set())
    print(f"{cid:<58}{len(fams):>11}   {sorted(fams)}")

n3 = sum(1 for cid in candidate_ids if len(marker_families.get(cid, set())) >= 3)
print(f"\nCandidates with >=3 distinct marker genes: {n3} / {len(candidate_ids)}")

node8 = [c for c in candidate_ids if "NODE_8_" in c]
if node8:
    fams = marker_families.get(node8[0], set())
    print(f"\nNODE_8 specifically: {len(fams)} distinct markers -> {sorted(fams)}")
