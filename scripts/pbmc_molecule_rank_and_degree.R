library(Seurat)
library(here)
library(pixelatorR)
library(ggplot2)
library(dplyr)
library(tidygraph)

here::i_am("scripts/pbmc_molecule_rank_and_degree.R")

pg_data <- ReadPNA_Seurat(here(
    "pbmc/data_untracked/PNA062_unstim_PBMCs_1000cells_S02_S2.layout.pxl"
))

molecule_rank_plot <- MoleculeRankPlot(pg_data) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    labs(title = "PBMC")
ggsave(
    here("pbmc/output/molecule_rank_plot.pdf"),
    molecule_rank_plot,
    width = 4,
    height = 4
)

## Degree distribution plot
set.seed(112)
sample_cells <- sample(colnames(pg_data), size = 100)
pg_data_sampled <- subset(pg_data, cells = sample_cells)

# Load cell graphs
pg_data_sampled <- LoadCellGraphs(pg_data_sampled)
cg_list <- CellGraphs(pg_data_sampled)

# Compute linker scores
calculate_degree <- function(cg) {
    cg@cellgraph |>
        activate(nodes) |>
        mutate(degree = centrality_degree()) |>
        as_tibble()
}

degree_list <- lapply(cg_list, calculate_degree)
degree_df <- bind_rows(degree_list, .id = "sample")

median_degree <- median(degree_df$degree)

degree_distribution <- degree_df |>
    group_by(degree) |>
    summarise(count = n()) |>
    mutate(percentage = count / sum(count) * 100) |>
    ggplot(aes(x = degree, y = percentage)) +
    geom_bar(stat = "identity", fill = "#69b3a2", alpha = 0.7) +
    geom_vline(xintercept = median_degree, linetype = "dashed") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    labs(
        x = "Degree",
        y = "Percentage of total (%)",
        title = "Degree Distribution (100 PBMC cells)"
    ) +
    annotate(
        "text",
        x = median_degree,
        y = Inf,
        label = sprintf("Median = %.2f", median_degree),
        vjust = 2,
        hjust = -0.1,
        size = 5,
        color = "black"
    )
degree_distribution

ggsave(
    here("pbmc/output/degree_distribution.pdf"),
    degree_distribution,
    width = 4,
    height = 4
)
