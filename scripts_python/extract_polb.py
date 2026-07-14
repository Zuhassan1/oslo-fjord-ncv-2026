import re, os
from collections import defaultdict

BEREN_ROOT = "/hdd/zoh/oslo_virus_project/giant_virus_screening"
PASSING_IDS = "passing_contig_ids.txt"
TARGET_FAMILY = "GVOGm0054"  # confirmed = PolB

with open(PASSING_IDS) as f:
    candidate_ids = [l.strip() for l in f if l.strip()]

by_sample = defaultdict(dict)
for cid in candidate_ids:
    m = re.match(r'^(ERR\d+)_(.+)$', cid)
    sample, node = m.group(1), m.group(2)
    by_sample[sample][node] = cid

def parse_fasta(path):
    header, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)

found = {}  # cid -> (seq, header) keeping the longest if multiple copies

for sample, node_to_cid in by_sample.items():
    marker_faa = f"{BEREN_ROOT}/{sample}_BEREN/Final_Results/NCLDV_Markers.faa"
    if not os.path.exists(marker_faa):
        continue
    for header, seq in parse_fasta(marker_faa):
        parts = header.split()
        if len(parts) < 2:
            continue
        marker_field, contig_field = parts[0], parts[1]
        if not marker_field.startswith(f"prodigal_{TARGET_FAMILY}."):
            continue
        node_id = re.sub(r'_\d+$', '', contig_field)
        if node_id in node_to_cid:
            cid = node_to_cid[node_id]
            if cid not in found or len(seq) > len(found[cid][0]):
                found[cid] = (seq, header)

print(f"Candidates with PolB (GVOGm0054): {len(found)} / {len(candidate_ids)} (expect 5)")
with open("polb_candidates.faa", "w") as out:
    for cid, (seq, orig_header) in found.items():
        out.write(f">{cid}\n{seq}\n")
        print(f"  {cid}  (length {len(seq)} aa)")
