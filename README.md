# Oslo Fjord Nucleocytoviricota (NCV) Discovery — Reproducibility Package

This repository contains the custom analysis scripts, intermediate data tables, and final figures underlying the manuscript:

**"Genome-Resolved Discovery of Nucleocytoviricota in an Urban Nordic Fjord: A Short-Read Metagenomic Survey Reveals a Putatively Complete Asfuvirales-Related Genome and Evidence of Site-Associated Viral Recurrence"**

Full methodological details, including every disclosed parameter and tool version, are provided in the manuscript's Materials and Methods section.

> **Note on this release:** this version (v1.1.0) incorporates a set of corrections identified during an independent post-hoc audit of the manuscript against this repository's own underlying data, conducted after the v1.0.0 release. Every correction was independently re-verified against primary source files (raw tool outputs, run logs) before being applied — see "What changed in v1.1.0" below for the complete list. Consistent with this project's approach throughout, corrections are disclosed openly rather than silently folded in.

## What this repository is (and is not)

This repository contains the **custom code written for this study** — the scripts used to filter, tally, extract, and visualize results from third-party tool outputs. It does **not** contain:
- Raw sequencing reads (publicly available at ENA, accessions ERR14872774–ERR14872782, study ERP171883)
- Third-party tool source code (BEREN, TIGTOG, GVClass, NuPhylo, dRep, eggNOG-mapper, ARAGORN, etc. — see manuscript References for citations and original repositories)
- Large reference databases (Kraken2 PlusPFP, NCBI nr, eggNOG-DB) used during screening and annotation

## Repository structure

```
scripts_python/          Custom Python scripts (contamination screening, completeness
                          assessment, PolB extraction, tRNA/COG tallying, abundance
                          matrix construction, master table consolidation)
    superseded/           Earlier, imperfect versions of five scripts, retained
                          openly rather than deleted -- see "On the superseded
                          folders" below.
scripts_R/                Six R scripts (01-06), each independently reproducing
                          one main-text figure from the data tables in data_tables/
    superseded/           The original version of script 06, superseded after a
                          confirmed taxonomy-parsing bug -- see below.
figures/                  Final figures (PDF + PNG), as they appear in the manuscript
data_tables/              Intermediate data files required to run the R scripts
                          without re-running the full upstream pipeline, including
                          the master supplementary table (all 30 pre-filter
                          candidates' full screening/completeness/taxonomy record)
```

## On the superseded folders

**`scripts_python/superseded/`** — four scripts representing earlier versions of analyses found, during the study's own internal quality-control process, to contain errors:
- `summarize_abundance.py` and `summarize_abundance_highconf.py` — an earlier abundance-summary approach later found to be confounded by uneven candidate-catalog composition across origin samples; superseded by `summarize_abundance_final.py`.
- `tally_cog_categories.py` — run on an incorrectly-scoped 30-candidate protein set (rather than the final, validated 26); superseded by `tally_cog_categories_VALID26.py` in the main `scripts_python/` folder.
- `tally_trna.py` — an earlier tRNA/intron parser that did not correctly detect intron annotations embedded inline within a gene-entry line; superseded by `tally_trna_v2.py`.

**`scripts_R/superseded/`** — one script, added in v1.1.0:
- `06_taxonomy_agreement.R` (original) — extracted the first `NCLDV__`-prefixed substring from GVClass's raw, unweighted percentage breakdown to determine each candidate's taxonomic order, without checking whether that substring actually represented GVClass's own confident consensus call. This incorrectly credited GVClass with an order-level call it had not actually made for two candidates (including the study's complete genome, ERR14872777_NODE_8), inflating the reported TIGTOG–GVClass agreement count from the correct 14/26 to an incorrect 15/26. Superseded by the current `06_taxonomy_agreement.R`, which reads GVClass's own `taxonomy_majority` consensus field directly.

These are retained, rather than deleted, as an open record of the corrections made during analysis. **Only the scripts in the top level of each `scripts_*/` folder (not `superseded/`) were used to produce the results and figures reported in the current manuscript version.**

## Reproducing the figures

**Requirements:** R (tested on R 4.6.1) with packages: `ape`, `tidytree`, `treeio`, `ggtree`, `ggplot2`, `dplyr`, `pheatmap`, `stringr`.

**Important — run each script in its own fresh R session.** During development, running multiple `ggtree`/`tidytree`-dependent scripts sequentially within a single persistent R session (e.g., pasting one script after another into an interactive RStudio console) was found to occasionally trigger a package-namespace conflict unrelated to script correctness. Each script below has been independently verified to run correctly as a standalone process. To reproduce:

```bash
cd data_tables/
Rscript ../scripts_R/01_polb_tree_main.R
Rscript ../scripts_R/02_polb_node8_inset.R
Rscript ../scripts_R/03_abundance_heatmap.R
Rscript ../scripts_R/04_huk_persistence_panel.R
Rscript ../scripts_R/05_cog_categories.R
Rscript ../scripts_R/06_taxonomy_agreement.R
```

Run each command as a **separate** invocation (each `Rscript` call starts a genuinely fresh R process) rather than sourcing all six inside one interactive session. Each script reads its required input file(s) directly from the current working directory (`data_tables/`) and writes its output figure(s) to the same location.

Script 5 (`05_cog_categories.R`) independently re-derives COG category counts directly from the raw eggNOG-mapper annotation file, rather than using precomputed values, as a genuine reproducibility check; its output was confirmed to match the manuscript's reported counts exactly (160/69/54/39/29/19/16/15/14/12/11/11/9/6/5/5/5/4/2/2/2 for categories S/L/K/O/F/A/G/M/J/H/T/E/I/C/B/P/V/U/D/Q/W respectively).

Script 6 (`06_taxonomy_agreement.R`) reads GVClass's own `taxonomy_majority` consensus field to determine each candidate's order-level call, rather than parsing the raw percentage breakdown, correctly identifying two candidates (including ERR14872777_NODE_8) for which GVClass did not reach a confident order-level consensus. This reproduces the manuscript's corrected agreement count exactly (14 Agree / 2 Both: no signal / 2 GVClass: no confident order call / 8 Disagree).

| Script | Figure produced | Required input file(s) |
|---|---|---|
| `01_polb_tree_main.R` | Figure 1 | `allseqs.nwk` |
| `02_polb_node8_inset.R` | Figure 2 | `allseqs.nwk` |
| `03_abundance_heatmap.R` | Figure 3 | `abundance_matrix_RPM.tsv` |
| `04_huk_persistence_panel.R` | Figure 4 | `abundance_matrix_RPM.tsv` |
| `05_cog_categories.R` | Figure 5 | `candidate_annotations_VALID26.emapper.annotations` |
| `06_taxonomy_agreement.R` | Figure 6 + supplementary table | `tigtog_results_v5.prediction_result.tsv`, `gvclass_summary.tsv` |

## Python scripts

The Python scripts in `scripts_python/` were run against outputs of the third-party tools listed in the manuscript's Methods section (BEREN, Kraken2, DIAMOND, dRep, eggNOG-mapper, ARAGORN), at the pipeline stages described there. They are provided for transparency and to allow inspection of the exact filtering/parsing logic used at each step, rather than as a fully automated one-command pipeline (running them end-to-end requires the upstream tool outputs, which are not included in this repository due to size, but are reproducible from the raw reads using the tool versions and parameters disclosed in the manuscript).

`build_master_table.py` (new in v1.1.0) consolidates the complete per-candidate decision chain for all 30 pre-filter candidates — contamination screening (Kraken2/DIAMOND percentages), ViralRecall score, rRNA screening, marker gene count and identity, circularity, terminal inverted repeat status, TIGTOG and GVClass taxonomic calls, and final abundance-analysis inclusion status — into `data_tables/candidate_master_supplementary_table.tsv`. It reuses the exact, independently-verified logic of `count_markers.py`, `check_circularity.py`, `parse_tir.py`, and `compute_contamination.py`, restructured to produce one consolidated table rather than four separate printed outputs.

## What changed in v1.1.0

Following an independent post-hoc audit comparing the manuscript against this repository's own data, the following corrections were made, each verified against primary source files before being applied:

1. **NODE_8 taxonomic framing** (Abstract, Results, Discussion, Figure 6): corrected to accurately reflect that GVClass did not reach a confident order-level consensus for this candidate (raw evidence: 25% Asfuvirales, 25% uninformative, 16.7% Pimascovirales at order level), rather than stating both classifiers agreed. TIGTOG's independent call (Asfuvirales, 0.65 confidence) and PolB phylogenetic placement remain the basis for this candidate's taxonomic discussion.
2. **Taxonomic agreement statistic**: corrected from 15/26 (57.7%) to 14/26 (53.8%), using GVClass's own `taxonomy_majority` field rather than its raw percentage breakdown.
3. **Candidate catalog length**: corrected from 1,643,672 bp to 1,485,520 bp, verified independently via `seqkit stats` and direct summation of the FASTA file.
4. **Huk home-dominance ratio range**: corrected from 0.9–1.4× to 0.93–1.11×, recomputed directly from the abundance matrix for the four ANI-confirmed candidates specifically.
5. **Dereplication algorithm**: corrected from an assumed "whole-genome ANI (ANImf)" to the actual algorithm used, fastANI (confirmed from the dRep run log), with real bidirectional alignment-coverage values (92–100%) added.
6. **Flood-event description**: the "<12 hours" timing claim was confirmed directly against the companion study's own text and retained; the more specific "combined sewer overflow" mechanism, which the companion study does not state, was removed.
7. Wording precision improvements: "recently published" → "preprint" (the reference methodology remains a bioRxiv preprint); HGT speculation on excluded candidates hedged as one possible explanation rather than asserted; "recurring exclusively" clarified to mean independently-assembled genomes (reads from these populations are also detected, at lower abundance, at Operastranda); "complete genome" → "putatively complete" at the two primary claim points; SNOWGIANTS grant number corrected (101150901 → 101210192, verified via CORDIS).
8. Title updated to match the corrected completeness and persistence framing.
9. Master supplementary table (`candidate_master_supplementary_table.tsv`) added, consolidating the full 30-candidate decision chain that was previously only reconstructable by running multiple separate scripts.

A sensitivity check confirmed the flood-candidate home-dominance ratio range (6.9–8.5×) is unchanged whether or not ERR14872777_NODE_333 (a candidate with conflicting taxonomic signal between classifiers) is included.

## Data availability

- Raw sequencing reads: ENA accessions ERR14872774–ERR14872782 (study ERP171883)
- Candidate genome assemblies: [ENA/NCBI accession to be added upon deposition]
- This repository: archived at Zenodo, DOI [10.5281/zenodo.21353320](https://doi.org/10.5281/zenodo.21353320) (concept DOI; always resolves to the latest version). This release (v1.1.0) is separately and permanently archived at [DOI to be added after this release is published — see instructions provided alongside this file].

## Citation

If you use these scripts or data, please cite the manuscript above and this repository's Zenodo DOI (10.5281/zenodo.21353320).

## Contact

[Author contact to be added]
