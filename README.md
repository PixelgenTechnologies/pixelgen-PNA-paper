# pixelgen-PNA-paper

Analysis code accompanying the Proximity Network Assay (PNA) manuscript *Cells as Single-Molecule Protein Proximity Networks*. PNA builds single-cell, single-molecule spatial networks of cell-surface proteins from DNA-barcoded antibodies and sequencing, enabling abundance, clustering, and colocalization analyses without optical imaging.

Preprint: [https://www.biorxiv.org/content/10.1101/2025.06.19.660329v3](https://www.biorxiv.org/content/10.1101/2025.06.19.660329v3)

## Repository structure

- `data/` – sample metadata and annotation tables used by the notebooks
- `scripts/` – main analysis notebooks (PBMC, SLE, PHA/stimulation, shared helpers)
- `raji/` – Raji cell-line analyses and MPX comparison scripts
- `car_t/` – CAR-T / Raji coculture analyses
- `roche/` – Illumina vs nanopore (SBX) comparison notebook
- `results/` – local outputs 

## Requirements

Install packages in your own R environment.

- [pixelatorR](https://github.com/PixelgenTechnologies/pixelatorR) v0.18.3
- Seurat v5.2.1
- tidyverse and related plotting packages (`ggplot2`, `ComplexHeatmap`, `patchwork`, etc.)

Primary sequencing data were processed with [nf-core/pixelator](https://github.com/nf-core/pixelator) to produce `.pxl` files.

## Data

Public datasets for this work are available at:

[https://software.pixelgen.com/datasets/pna-methods-paper](https://software.pixelgen.com/datasets/pna-methods-paper)

Download FASTQ and/or PXL files from that page and place them where the notebooks expect them (typically under `data/` or experiment-specific folders such as `raji/data_untracked/`). Sample identifiers and filenames are listed in [`data/sample_metadata.csv`](data/sample_metadata.csv).

## Getting started

```bash
git clone https://github.com/PixelgenTechnologies/pixelgen-PNA-paper.git
cd pixelgen-PNA-paper
```

Then open the repository in R / RStudio / Quarto and install the packages above.

### Suggested entry points

These notebooks are the most self-contained relative paths for public PXLs once files are in place:

| Analysis | Script |
| --- | --- |
| PBMC annotation (Figure 3) | [`scripts/pbmc_analysis.qmd`](scripts/pbmc_analysis.qmd) |
| Raji proximity scores (Figure 2) | [`raji/scripts/analysis.qmd`](raji/scripts/analysis.qmd) |
| PHA / resting PBMC | [`scripts/stim_analysis.qmd`](scripts/stim_analysis.qmd) |

Additional notebooks (CAR-T, combined paper script, Roche SBX, SLE, MPX comparison) live under `car_t/`, `scripts/`, `roche/`, and `raji/scripts/`. 

## Citation

If you use this code or the associated datasets, please cite:

> Cells as Single-Molecule Protein Proximity Networks  
> Filip Karlsson, Michele Simonetti, Christina Galonska, Max Jonatan Karlsson, Hanna van Ooijen, Divya Thiagarajan, Tomasz Kallas, Maud Schweitzer, Ludvig Larsson, Vincent van Hoef, Pouria Tajvar, Johan Dahlberg, Florian De Temmerman, Louise Leijonancker, Vanessa Trombin, Sylvain Geny, Rikard Forlin, Erika Negrini, Beijing Wu, Liu Xi, Stefan Petkov, Lovisa Franzen, Jessica Bunz, Christine Moge, Henrik Everberg, Petter Brodin, Alvaro Martinez Barrio, Simon Fredriksson  
> bioRxiv 2025.06.19.660329; doi: [https://doi.org/10.1101/2025.06.19.660329](https://doi.org/10.1101/2025.06.19.660329)

## License

This project is licensed under the [GPL-2.0 License](LICENSE).

## Contact

For questions about this repository, open a GitHub issue.
