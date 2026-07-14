import re

def get_candidates_in_faa(path):
    candidates = set()
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            seqid = line[1:].split()[0]          # true ID only, ignore # coord metadata
            candidate = re.sub(r'_\d+$', '', seqid)  # strip trailing _<gene_number>
            candidates.add(candidate)
    return candidates

with open("passing_contig_ids.txt") as f:
    valid_26 = set(l.strip() for l in f if l.strip())

in_faa = get_candidates_in_faa("candidate_proteins.faa")

print(f"Candidates represented in candidate_proteins.faa: {len(in_faa)}")
print(f"Valid, final 26-candidate list: {len(valid_26)}")

invalid = in_faa - valid_26
print(f"\nCandidates in the protein file that are NOT in the valid 26: {len(invalid)}")
for c in sorted(invalid):
    print(f"  {c}")

missing = valid_26 - in_faa
print(f"\nValid candidates missing from the protein file (should be none): {len(missing)}")
for c in sorted(missing):
    print(f"  {c}")
