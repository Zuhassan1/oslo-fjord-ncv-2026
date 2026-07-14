import re, sys
from collections import defaultdict

NODES_DMP = "/hdd/zoh/databases/diamond_nr/nodes.dmp"
DIAMOND_OUT = "diamond_screening_results.tsv"
MARKER_FAA = "/hdd/zoh/oslo_virus_project/giant_virus_screening/ERR14872777_BEREN/Final_Results/NCLDV_Markers.faa"

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

REGION1_END = 127019
REGION2_START = 128023

# Parse protein coordinates
proteins = {}
with open("node8_protein_coords.txt") as f:
    for line in f:
        m = re.match(r'^>(\S+) # (\d+) # (\d+) #', line)
        if m:
            pid, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            proteins[pid] = (start, end)

def region_of(start, end):
    mid = (start + end) // 2
    if mid <= REGION1_END:
        return "Region1"
    elif mid >= REGION2_START:
        return "Region2"
    return "Gap"

region_counts = defaultdict(lambda: {"n_prot": 0, "n_hits": 0, "n_cellular": 0, "n_viral": 0})
for pid, (s, e) in proteins.items():
    region_counts[region_of(s, e)]["n_prot"] += 1

# DIAMOND hits for NODE_8 proteins
with open(DIAMOND_OUT) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        qseqid, staxid_field = fields[0], fields[6]
        if qseqid not in proteins:
            continue
        s, e = proteins[qseqid]
        r = region_of(s, e)
        region_counts[r]["n_hits"] += 1
        cls = classify(staxid_field.split(";")[0].strip())
        if cls == "cellular":
            region_counts[r]["n_cellular"] += 1
        elif cls == "viral":
            region_counts[r]["n_viral"] += 1

# Marker gene positions
marker_hits = []
with open(MARKER_FAA) as f:
    for line in f:
        if not line.startswith(">"):
            continue
        parts = line[1:].split()
        if len(parts) < 2 or "NODE_8_length_176732" not in parts[1]:
            continue
        marker_field, contig_field = parts[0], parts[1]
        m2 = re.match(r'^prodigal_(GVOGm\d+)\.copy\d+$', marker_field)
        if not m2:
            continue
        pid = contig_field.split(" # ")[0] if " # " in contig_field else None
        gene_num = re.search(r'_(\d+)$', contig_field.split()[0] if not pid else pid)

print("=== Region breakdown ===")
for r in ["Region1", "Region2", "Gap"]:
    c = region_counts[r]
    pct = (c["n_cellular"] / c["n_hits"] * 100) if c["n_hits"] else 0
    print(f"{r}: {c['n_prot']} proteins, {c['n_hits']} DIAMOND hits, {c['n_cellular']} cellular ({pct:.1f}%), {c['n_viral']} viral")

print("\n=== Marker gene protein IDs on NODE_8 (raw, for manual position lookup) ===")
with open(MARKER_FAA) as f:
    for line in f:
        if line.startswith(">") and "NODE_8_length_176732" in line:
            print(line.strip())
