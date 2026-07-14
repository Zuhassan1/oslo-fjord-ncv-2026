from collections import defaultdict

ANNOT_FILE = "candidate_annotations_VALID26.emapper.annotations"

cog_counts = defaultdict(int)
total_annotated = 0
total_lines = 0

with open(ANNOT_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        total_lines += 1
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue
        cog_field = fields[6].strip()
        if cog_field and cog_field != "-":
            total_annotated += 1
            for letter in cog_field:
                cog_counts[letter] += 1

print(f"Total protein lines in annotation file: {total_lines}")
print(f"Proteins with a COG category assigned: {total_annotated}")
print(f"\n{'COG':<6}{'Count':>8}")
for letter, count in sorted(cog_counts.items(), key=lambda x: -x[1]):
    print(f"{letter:<6}{count:>8}")
print(f"\nTotal COG letter-assignments: {sum(cog_counts.values())}")
