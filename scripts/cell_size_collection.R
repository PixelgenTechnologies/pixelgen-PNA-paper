library(pixelatorR)
library(here)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(purrr)
library(stringr)
library(forcats)
library(readr)
library(readxl)
here::i_am("scripts/cell_size_collection.R")


raji_data <- ReadPNA_Seurat(here(
    "raji/data_untracked/PNA062_Raji_1000cells_S06_S6.layout.pxl"
))


mrp_raji <- MoleculeRankPlot(raji_data) +
    theme_bw() +
    labs(title = "Raji Cells")
ggsave(here("raji/output/mrp_raji.pdf"), mrp_raji, width = 8, height = 5)

raji_cell_sizes <- raji_data[[]] %>%
    summarise(cell_size = mean(n_umi), sd = sd(n_umi)) |>
    mutate(dataset = "raji")
