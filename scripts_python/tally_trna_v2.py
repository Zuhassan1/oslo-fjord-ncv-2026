import re

with open("/hdd/zoh/oslo_virus_project/giant_virus_screening/candidate_catalog/passing_contig_ids.txt") as f:
    valid_26 = [l.strip() for l in f if l.strip()]

gene_counts = {}
trna_types = {}
intron_details = {}
current = None

with open("/hdd/zoh/oslo_virus_project/giant_virus_screening/trna_scan/aragorn_results.txt") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">") and not line.startswith(">end"):
            current = line[1:].strip()
            trna_types[current] = []
            intron_details[current] = []
        elif "genes found" in line or "gene found" in line:
            n = int(line.split()[0])
            gene_counts[current] = n
        elif re.match(r'^\s*\d+\s+tRNA-', line):
            species = line.split()[1].replace("tRNA-", "")
            trna_types[current].append(species)
            # Check for an inline intron annotation anywhere on this same line
            m = re.search(r'i\((\d+),(\d+)\)', line)
            if m:
                intron_details[current].append((species, int(m.group(1)), int(m.group(2))))

print(f"Candidates found: {len(gene_counts)} (expect 26)")
print(f"\n{'Candidate':<58}{'tRNAs':>7}{'Introns':>9}   Types")
total_trna, total_introns = 0, 0
for c in valid_26:
    n = gene_counts.get(c, 0)
    introns = intron_details.get(c, [])
    total_trna += n
    total_introns += len(introns)
    print(f"{c:<58}{n:>7}{len(introns):>9}   {sorted(set(trna_types.get(c, [])))}")

print(f"\nTotal tRNA genes: {total_trna}")
print(f"Total introns detected: {total_introns}")
print(f"Candidates with >=1 tRNA: {sum(1 for c in valid_26 if gene_counts.get(c,0) > 0)} / 26")
print(f"\nIntron details:")
for c, introns in intron_details.items():
    for species, pos, length in introns:
        print(f"  {c}: tRNA-{species}, intron at position {pos}, length {length} bp")
