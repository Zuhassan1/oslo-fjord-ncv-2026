# Oslo Fjord Nucleocytoviricota (NCV) Discovery — Reproducibility Package

This repository contains the custom analysis scripts, intermediate data tables, and final figures underlying the manuscript:

**"Genome-Resolved Discovery of Nucleocytoviricota in an Urban Nordic Fjord: A Short-Read Metagenomic Survey Reveals a Complete Asfuvirales-Related Genome and Evidence of a Site-Persistent Viral Population"**

Full methodological details, including every disclosed parameter and tool version, are provided in the manuscript's Materials and Methods section and in the accompanying `materials_and_methods_validated.md` supplementary record.

## What this repository is (and is not)

This repository contains the **custom code written for this study** — the scripts used to filter, tally, extract, and visualize results from third-party tool outputs. It does **not** contain:
- Raw sequencing reads (publicly available at ENA, accessions ERR14872774–ERR14872782, study ERP171883)
- Third-party tool source code (BEREN, TIGTOG, GVClass, NuPhylo, dRep, eggNOG-mapper, ARAGORN, etc. — see manuscript References for citations and original repositories)
- Large reference databases (Kraken2 PlusPFP, NCBI nr, eggNOG-DB) used during screening and annotation

## Repository structure

```
scripts_python/          Custom Python scripts (contamination screening, completeness
                          assessment, PolB extraction, tRNA/COG tallying, abundance
                          matrix construction)
    superseded/           Earlier, imperfect versions of four scripts, retained
                          openly rather than deleted -- see "On the superseded/
                          folder" below.
scripts_R/                Six R scripts (01-06), each independently reproducing
                          one main-text figure from the data tables in data_tables/
figures/                  Final figures (PDF + PNG), as they appear in the manuscript
data_tables/               Intermediate data files required to run the R scripts
                          without re-running the full upstream pipeline
```

## On the superseded/ folder

Four scripts in `scripts_python/superseded/` represent earlier versions of analyses that were found, during the study's own internal quality-control process, to contain errors:
- `summarize_abundance.py` and `summarize_abundance_highconf.py` — an earlier abundance-summary approach later found to be confounded by uneven candidate-catalog composition across origin samples; superseded by `summarize_abundance_final.py`.
- `tally_cog_categories.py` — run on an incorrectly-scoped 30-candidate protein set (rather than the final, validated 26); superseded by `tally_cog_categories_VALID26.py` in the main `scripts_python/` folder.
- `tally_trna.py` — an earlier tRNA/intron parser that did not correctly detect intron annotations embedded inline within a gene-entry line; superseded by `tally_trna_v2.py`.

These are retained, rather than deleted, as an open record of the corrections made during analysis, consistent with the manuscript's approach of explicitly disclosing methodological deviations and corrections rather than silently resolving them. **Only the scripts in the top level of `scripts_python/` (not `superseded/`) were used to produce the results and figures reported in the manuscript.**

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

## Data availability

- Raw sequencing reads: ENA accessions ERR14872774–ERR14872782 (study ERP171883)
- Candidate genome assemblies: [ENA/NCBI accession to be added upon deposition]
- This repository: archived at Zenodo under DOI [to be added upon first release]

## Citation

If you use these scripts or data, please cite the manuscript above and this repository's Zenodo DOI.

## Contact

[Author contact to be added]
