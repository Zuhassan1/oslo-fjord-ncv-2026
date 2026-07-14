import re

with open("/hdd/zoh/oslo_virus_project/giant_virus_screening/candidate_catalog/passing_contig_ids.txt") as f:
    valid_26 = [l.strip() for l in f if l.strip()]

gene_counts = {}
trna_types = {}
intron_counts = {}
current = None

with open("/hdd/zoh/oslo_virus_project/giant_virus_screening/trna_scan/aragorn_results.txt") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            current = line[1:].strip()
            trna_types[current] = []
            intron_counts[current] = 0
        elif "genes found" in line or "gene found" in line:
            n = int(line.split()[0])
            gene_counts[current] = n
        elif re.match(r'^\s*\d+\s+tRNA-', line):
            species = line.split()[1].replace("tRNA-", "")
            trna_types[current].append(species)
        elif line.strip().startswith("i("):
            intron_counts[current] += 1

print(f"Candidates found in ARAGORN output: {len(gene_counts)} (expect 26)")
missing = set(valid_26) - set(gene_counts.keys())
print(f"Missing from output (should be none): {missing}\n")

print(f"{'Candidate':<58}{'tRNAs':>7}{'Introns':>9}   Types")
total_trna = 0
total_introns = 0
for c in valid_26:
    n = gene_counts.get(c, 0)
    intr = intron_counts.get(c, 0)
    types = trna_types.get(c, [])
    total_trna += n
    total_introns += intr
    print(f"{c:<58}{n:>7}{intr:>9}   {sorted(set(types))}")

print(f"\nTotal tRNA genes across all 26 candidates: {total_trna}")
print(f"Total introns detected: {total_introns}")
print(f"Candidates with >=1 tRNA: {sum(1 for c in valid_26 if gene_counts.get(c,0) > 0)} / 26")
