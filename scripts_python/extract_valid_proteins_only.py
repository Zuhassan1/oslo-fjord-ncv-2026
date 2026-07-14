import re

def parse_fasta(path):
    header, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]  # keep full header (includes coordinate metadata)
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)

with open("passing_contig_ids.txt") as f:
    valid_26 = set(l.strip() for l in f if l.strip())

kept, dropped = 0, 0
with open("candidate_proteins_VALID26.faa", "w") as out:
    for header, seq in parse_fasta("candidate_proteins.faa"):
        seqid = header.split()[0]
        candidate = re.sub(r'_\d+$', '', seqid)
        if candidate in valid_26:
            out.write(f">{header}\n{seq}\n")
            kept += 1
        else:
            dropped += 1

print(f"Kept: {kept} proteins from valid candidates")
print(f"Dropped: {dropped} proteins from the 4 excluded candidates")
