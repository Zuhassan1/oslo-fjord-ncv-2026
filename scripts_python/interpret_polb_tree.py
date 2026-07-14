from Bio import Phylo
import io

TREE_FILE = "polb_tree_results/allseqs.nwk"
OUR_CANDIDATES = [
    "ERR14872774_NODE_186_length_52281_cov_9.344311",
    "ERR14872777_NODE_111_length_52561_cov_10.705310",
    "ERR14872777_NODE_118_length_50144_cov_20.226856",
    "ERR14872777_NODE_20_length_137217_cov_20.093298",
    "ERR14872777_NODE_8_length_176732_cov_16.217617",
]

tree = Phylo.read(TREE_FILE, "newick")
all_terminals = tree.get_terminals()
print(f"Total terminal tips in tree: {len(all_terminals)} (expect 113)\n")

for cand in OUR_CANDIDATES:
    target = None
    for term in all_terminals:
        if term.name == cand:
            target = term
            break
    if target is None:
        print(f"!!! {cand}: NOT FOUND IN TREE")
        continue

    path = tree.get_path(target)
    parent = path[-2] if len(path) >= 2 else tree.root
    siblings = [t.name for t in parent.get_terminals() if t.name != cand]

    print(f"=== {cand} ===")
    print(f"  Immediate clade contains {len(siblings)+1} tips total")
    print(f"  Nearest neighbors (same immediate clade): {siblings[:8]}")
    if hasattr(parent, "confidence") and parent.confidence is not None:
        print(f"  Support at this node: {parent.confidence}")
    print()
