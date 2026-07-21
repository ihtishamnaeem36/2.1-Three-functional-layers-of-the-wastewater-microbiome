# The acquired resistome of untreated wastewater is inversely linked to native antibiotic-biosynthesis potential

**A multi-source metagenomic study from Khyber Pakhtunkhwa, Pakistan**

[![DOI](https://zenodo.org/badge/1305511427.svg)](https://doi.org/10.5281/zenodo.21439016)

This repository holds the code, processed data, and figures for a shotgun-metagenomic study of wastewater from three cities in Khyber Pakhtunkhwa (Swat, Mardan, Peshawar), sampled across hospital, slaughterhouse, and community sites (n = 18). The study treats the community as three functional layers — metabolism, the intrinsic resistome, and the acquired resistome — and asks how each is organised and how they interact.


## Repository layout

```
github/
├── README.md                 this file
├── LICENSE                   MIT licence (code)
├── CITATION.cff              how to cite this work
├── requirements.txt          Python dependencies (pip)
├── environment.yml           Python dependencies (conda)
├── .gitignore
├── code/
│   └── analysis_pipeline.py  full analysis: metabolism, CARD resistome,
│                             KEGG intrinsic resistance, biosynthesis screen, mobilome
├── data/
│   ├── processed/            result tables (CSV/TSV): diversity, PERMANOVA,
│   │                         PCoA, network metrics, resistome and correlation tables
│   ├── mobilome/             per-sample mobile-element markers
│   └── README.md             data dictionary and availability
├── figures/
│   ├── main/                 Figures 1–9 (PNG at 600 dpi + SVG vector)
│   └── supplementary/        Supplementary Figures S1–S4
└── docs/
    ├── Methods.txt           full methods
    ├── Results.txt           full results narrative
    ├── Synopsis.md           project synopsis
    └── Figure_legends.txt    figure legends
```

## Data availability

Raw sequence reads are deposited in the NCBI Sequence Read Archive under **BioProject PRJNA1463247**, with the per-sample run accessions for the 18 samples listed in Supplementary Table S1 of the manuscript; the reads are held private until publication. This repository is archived on Zenodo: **https://doi.org/10.5281/zenodo.21439016**. Processed KO and ARG abundance tables and all result tables used to generate the figures are in `data/`. Large upstream inputs (the KEGG functional annotation table and the CARD alignment) are not committed here because of size; they are available on request or via the deposited BioProject.

## Reproducing the analysis

The pipeline is a single annotated script, `code/analysis_pipeline.py`, organised in five parts:

1. Metabolic pipeline — KO abundance, alpha/beta diversity, co-occurrence network, centralization.
2. CARD acquired resistome — ARG profiling, composition, burden, resistome–architecture tests.
3. KEGG intrinsic resistance — β-lactam, vancomycin, and cationic-antimicrobial-peptide resistance pathways.
4. Biosynthesis screen — variable functional categories, producer–acquirer axis.
5. Mobilome — class 1 integron and transposase markers, integration with the acquired resistome.

Each part reads its own input table; set the input and output paths at the top of each part before running. The script installs its own dependencies on first run, or you can create the environment beforehand:

```bash
# option A: pip
pip install -r requirements.txt

# option B: conda
conda env create -f environment.yml
conda activate wastewater-layers

python code/analysis_pipeline.py
```

Statistics use SciPy and scikit-bio conventions: PERMANOVA (999 permutations), Kruskal–Wallis, Mann–Whitney U, Spearman correlation, and Benjamini–Hochberg false-discovery correction. Networks use NetworkX (Spearman |r| ≥ 0.6, P < 0.001).

## Scope and limits

This is a first regional characterisation with a small, balanced design (n = 18, two replicates per city–environment cell). Within-cell power is low, the intrinsic-resistance PERMANOVA is borderline (P = 0.05), and biosynthesis potential is a gene-abundance proxy rather than measured production. Metagenomics reads potential, not expression, and short reads leave the integrons unassigned to host taxa. Results should be read in that light.

## Citation

If you use this code or data, please cite the work as described in `CITATION.cff`.

## Licence

Code is released under the MIT Licence (`LICENSE`). Processed data and figures are shared for reuse with attribution; please cite the associated publication.

## Contact

Author: Ihtisham Naeem (ihtishamnaeem36@gmail.com)

Corresponding authors: Ishaq Khan (ishaq@uswat.edu.pk); Muhammad Shafiq (drshafiq@stu.edu.cn).
