library(pixelatorR)
library(here)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)
library(ggrastr) # for geom_point_rast (rasterized points)
here::i_am("raji/scripts/score_stability.R")


########################################################
# Preprocess the PNA data
########################################################

# Read PNA data
pna_raji <- ReadPNA_Seurat(here(
    "raji/data_untracked/PNA062_Raji_1000cells_S06_S6.analysis.pxl"
))

# Sample 5 cells: one from each quintile of n_umi (0–20%, 20–40%, …, 80–100%)
set.seed(12)
sorted_cells <- names(sort(pna_raji$n_umi))
n <- length(sorted_cells)
# Split into 5 equal-sized bins and take the middle cell from each bin
bin_size <- n %/% 5
positions <- (seq_len(5) - 1) * bin_size + (bin_size %/% 2 + 1) # middle of each bin
positions <- pmin(positions, n) # clamp so we never exceed length
cell_selection <- sorted_cells[positions]
cell_selection_sizes <- pna_raji$n_umi[cell_selection]

# Density of n_umi as ggplot (for inspection / debugging)
# Indicate which cells are used in cell_selection on the plot
dens <- density(pna_raji$n_umi)
cell_selection_n_umi <- pna_raji$n_umi[cell_selection]

ggplot(data.frame(x = dens$x, y = dens$y), aes(x = x, y = y)) +
    geom_line() +
    geom_vline(
        xintercept = cell_selection_n_umi,
        linetype = "dashed",
        color = "red"
    ) +
    theme_bw() +
    labs(
        x = "n_umi",
        y = "Density",
        title = "Density of n_umi (red dashed: cell_selection)"
    ) +
    theme(legend.position = "none")


pna_raji_sampled <- subset(pna_raji, cells = cell_selection)

clust_scores <- ProximityScores(
    pna_raji_sampled,
    add_marker_counts = TRUE,
    add_marker_proportions = TRUE
)

pna_raji_sampled <- LoadCellGraphs(pna_raji_sampled)


# Repeat the proximity score computation 100 times with different seeds
# to assess stability; store each result (marker pairs + log2_ratio + component) in a list
res_50_list <- lapply(seq_len(100), function(i) {
    pixelatorRinternal::ComputeProximityScores(
        pna_raji_sampled,
        iterations = 50,
        seed = 123L + i, # different seed per run for independent replicates
        assay = "PNA",
        cells = cell_selection
    ) |>
        select(marker_1, marker_2, log2_ratio, component)
})
saveRDS(res_50_list, here("raji/data_untracked/res_50_list.rds"))
res_50_list <- readRDS(here("raji/data_untracked/res_50_list.rds"))

res_100_list <- lapply(seq_len(100), function(i) {
    pixelatorRinternal::ComputeProximityScores(
        pna_raji_sampled,
        iterations = 100,
        seed = 123L + i, # different seed per run for independent replicates
        assay = "PNA",
        cells = cell_selection
    ) |>
        select(marker_1, marker_2, log2_ratio, component)
})
saveRDS(res_100_list, here("raji/data_untracked/res_100_list.rds"))
res_100_list <- readRDS(here("raji/data_untracked/res_100_list.rds"))

res_500_list <- lapply(seq_len(100), function(i) {
    pixelatorRinternal::ComputeProximityScores(
        pna_raji_sampled,
        iterations = 500,
        seed = 123L + i,
        assay = "PNA",
        cells = cell_selection
    ) |>
        select(marker_1, marker_2, log2_ratio, component)
})
saveRDS(res_500_list, here("raji/data_untracked/res_500_list.rds"))
res_500_list <- readRDS(here("raji/data_untracked/res_500_list.rds"))

res_1000_list <- lapply(seq_len(100), function(i) {
    pixelatorRinternal::ComputeProximityScores(
        pna_raji_sampled,
        iterations = 1000,
        seed = 123L + i,
        assay = "PNA",
        cells = cell_selection
    ) |>
        select(marker_1, marker_2, log2_ratio, component)
})
saveRDS(res_1000_list, here("raji/data_untracked/res_1000_list.rds"))

res_250_list <- lapply(seq_len(100), function(i) {
    pixelatorRinternal::ComputeProximityScores(
        pna_raji_sampled,
        iterations = 250,
        seed = 123L + i,
        assay = "PNA",
        cells = cell_selection
    ) |>
        select(marker_1, marker_2, log2_ratio, component)
})
saveRDS(res_250_list, here("raji/data_untracked/res_250_list.rds"))

res_50_list <- readRDS(here("raji/data_untracked/res_50_list.rds"))
res_100_list <- readRDS(here("raji/data_untracked/res_100_list.rds"))
res_250_list <- readRDS(here("raji/data_untracked/res_250_list.rds"))
res_500_list <- readRDS(here("raji/data_untracked/res_500_list.rds"))
res_1000_list <- readRDS(here("raji/data_untracked/res_1000_list.rds"))


# Bind replicates and tag each dataset with its iteration count.
res_50_df <- bind_rows(res_50_list, .id = "rep") |>
    mutate(iterations = "50_permutations")
res_100_df <- bind_rows(res_100_list, .id = "rep") |>
    mutate(iterations = "100_permutations")
res_250_df <- bind_rows(res_250_list, .id = "rep") |>
    mutate(iterations = "250_permutations")
res_500_df <- bind_rows(res_500_list, .id = "rep") |>
    mutate(iterations = "500_permutations")
res_1000_df <- bind_rows(res_1000_list, .id = "rep") |>
    mutate(iterations = "1000_permutations")

# Bind all datasets together
res_all <- res_50_df |>
    bind_rows(res_100_df) |>
    bind_rows(res_250_df) |>
    bind_rows(res_500_df)

saveRDS(res_all, here("raji/data_untracked/res_all.rds"))
#res_all <- readRDS(here("raji/data_untracked/res_all.rds"))
########################################################
# Summarize stability: SD of log2_ratio across 100 replicates
# per (marker_1, marker_2, component, iterations)
########################################################
stability_summary <- res_all |>
    group_by(iterations, component, marker_1, marker_2) |>
    summarise(
        mean_log2_ratio = mean(log2_ratio, na.rm = TRUE),
        conf_95_lower = mean_log2_ratio - 1.96 * sd(log2_ratio, na.rm = TRUE),
        conf_95_upper = mean_log2_ratio + 1.96 * sd(log2_ratio, na.rm = TRUE),
        sd_log2_ratio = sd(log2_ratio, na.rm = TRUE),
        n_reps = n(),
        .groups = "drop"
    ) |>
    mutate(
        cell_size = cell_selection_sizes[match(component, cell_selection)]
    ) |>
    mutate(
        iterations = factor(
            iterations,
            levels = c(
                "50_permutations",
                "100_permutations",
                "250_permutations",
                "500_permutations"
            )
        )
    )


# Low abundant pair: select pair with mean counts in 100–200 range, present in all 5 cells
low_abundant_candidates <- clust_scores |>
    group_by(marker_1, marker_2) |>
    mutate(
        mean_count_1 = mean(count_1, na.rm = TRUE),
        mean_count_2 = mean(count_2, na.rm = TRUE)
    ) |>
    filter(abs(log2_ratio) > 0) |>
    filter(
        (mean_count_1 < 200 & mean_count_1 > 100) |
            (mean_count_2 < 200 & mean_count_2 > 100)
    ) |>
    filter(n() == 5) |>
    ungroup() |>
    distinct(marker_1, marker_2)
if (nrow(low_abundant_candidates) == 0) {
    stop("No low-abundant pair with n() == 5 found.")
}
# Prefer CD32/CD52 if present, otherwise use first candidate
cd32_cd52 <- low_abundant_candidates |>
    filter(marker_1 == "CD32" & marker_2 == "CD52")
if (nrow(cd32_cd52) > 0) {
    marker1_low <- "CD32"
    marker2_low <- "CD52"
} else {
    marker1_low <- low_abundant_candidates$marker_1[1]
    marker2_low <- low_abundant_candidates$marker_2[1]
}

# Add cell_1..cell_5 by cell size (smallest = cell_1, largest = cell_5) for the low-abundant pair
low_abundant_data <- stability_summary |>
    filter(marker_1 == marker1_low & marker_2 == marker2_low) |>
    mutate(
        cell_id = case_when(
            cell_size == min(cell_size) ~ "cell_1",
            cell_size == sort(unique(cell_size))[2] ~ "cell_2",
            cell_size == sort(unique(cell_size))[3] ~ "cell_3",
            cell_size == sort(unique(cell_size))[4] ~ "cell_4",
            cell_size == max(cell_size) ~ "cell_5",
            TRUE ~ as.character(cell_size)
        )
    ) |>
    mutate(
        cell_size_k = round(cell_size / 1000, 2),
        cell_size_label = paste0(cell_id, " (", cell_size_k, "k)")
    ) |>
    mutate(cell_id = factor(cell_id, levels = paste0("cell_", 1:5)))

# One plot per cell (cell_1 … cell_5), then concatenate into a single panel
plot_list_low <- lapply(paste0("cell_", 1:5), function(cid) {
    d <- low_abundant_data |> filter(cell_id == cid)
    if (nrow(d) == 0) {
        return(ggplot())
    }
    label <- d$cell_size_label[1]
    ggplot(d, aes(x = iterations, y = mean_log2_ratio, color = iterations)) +
        geom_point() +
        geom_errorbar(
            aes(
                ymin = mean_log2_ratio - sd_log2_ratio,
                ymax = mean_log2_ratio + sd_log2_ratio
            ),
            width = 0.15
        ) +
        theme_bw() +
        ylim(-1, 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
        labs(
            x = NULL,
            y = "Log2 Ratio",
            title = label
        ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            aspect.ratio = 1,
            plot.title = element_text(size = 12)
        )
})
p_low_marker <- wrap_plots(plot_list_low, nrow = 1) +
    plot_annotation(
        title = glue::glue(
            "Low-abundant pair {marker1_low}_{marker2_low}: stability per cell (cell_1 = smallest, cell_5 = largest)"
        ),
        theme = theme(plot.title = element_text(size = 14))
    )
ggsave(
    here(glue::glue(
        "raji/output/score_stability/{marker1_low}_{marker2_low}_by_cell.pdf"
    )),
    p_low_marker,
    width = 14,
    height = 4
)


marker1 <- "CD21"
marker2 <- "CD55"
p_medium_marker <- stability_summary |>
    filter(marker_1 == marker1 & marker_2 == marker2) |>
    mutate(
        cell_id = case_when(
            cell_size == min(cell_size) ~ "cell_1",
            cell_size == sort(unique(cell_size))[2] ~ "cell_2",
            cell_size == sort(unique(cell_size))[3] ~ "cell_3",
            cell_size == sort(unique(cell_size))[4] ~ "cell_4",
            cell_size == max(cell_size) ~ "cell_5",
            TRUE ~ as.character(cell_size)
        )
    ) |>
    mutate(
        cell_size_k = round(cell_size / 1000, 2), # convert to k, round to 2 decimals
        cell_size_label = paste0(cell_id, " (", cell_size_k, "k)")
    ) |>
    ggplot(aes(
        x = reorder(cell_size_label, cell_size),
        y = mean_log2_ratio,
        color = iterations,
        group = iterations
    )) +
    geom_point(position = position_dodge(width = 0.6)) +
    geom_errorbar(
        aes(
            ymin = mean_log2_ratio - sd_log2_ratio,
            ymax = mean_log2_ratio + sd_log2_ratio
        ),
        width = 0.2,
        position = position_dodge(width = 0.6)
    ) +
    theme_bw() +
    ylim(-1, 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(
        x = "Cell (1-5)",
        y = "Log2 Ratio",
        color = "Iterations",
        title = glue::glue(
            "Summary of 100 replicates of the proximity score computation for {marker1}_{marker2}"
        )
    ) +
    theme(
        aspect.ratio = 1,
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
    )
ggsave(
    here(glue::glue("raji/output/score_stability/{marker1}_{marker2}.pdf")),
    p_low_marker,
    width = 10,
    height = 5
)


stability_summary |>
    left_join(
        clust_scores |> select(marker_1, marker_2, component, count_1, count_2),
        by = c("marker_1", "marker_2", "component")
    ) |>
    filter(count_1 > 500 & count_2 > 500) |>
    filter(abs(mean_log2_ratio) > 0) |>
    group_by(marker_1, marker_2) |>
    filter(n() == 20) |>
    as.data.frame()
marker1 <- "CD19"
marker2 <- "CD19"
p_high_marker <- stability_summary |>
    filter(marker_1 == marker1 & marker_2 == marker2) |>
    mutate(
        cell_id = case_when(
            cell_size == min(cell_size) ~ "cell_1",
            cell_size == sort(unique(cell_size))[2] ~ "cell_2",
            cell_size == sort(unique(cell_size))[3] ~ "cell_3",
            cell_size == sort(unique(cell_size))[4] ~ "cell_4",
            cell_size == max(cell_size) ~ "cell_5",
            TRUE ~ as.character(cell_size)
        )
    ) |>
    mutate(
        cell_size_k = round(cell_size / 1000, 2), # convert to k, round to 2 decimals
        cell_size_label = paste0(cell_id, " (", cell_size_k, "k)")
    ) |>
    ggplot(aes(
        x = reorder(cell_size_label, cell_size),
        y = mean_log2_ratio,
        color = iterations,
        group = iterations
    )) +
    geom_point(position = position_dodge(width = 0.6)) +
    geom_errorbar(
        aes(
            ymin = mean_log2_ratio - sd_log2_ratio,
            ymax = mean_log2_ratio + sd_log2_ratio
        ),
        width = 0.2,
        position = position_dodge(width = 0.6)
    ) +
    theme_bw() +
    ylim(-1, 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(
        x = "Cell (1-5)",
        y = "Log2 Ratio",
        color = "Iterations",
        title = glue::glue(
            "Summary of 100 replicates of the proximity score computation for {marker1}_{marker2}"
        )
    ) +
    theme(
        aspect.ratio = 1,
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
    )
ggsave(
    here(glue::glue("raji/output/score_stability/{marker1}_{marker2}.pdf")),
    p_low_marker,
    width = 10,
    height = 5
)

p_all <- p_low_marker +
    p_medium_marker +
    p_high_marker +
    plot_layout(ncol = 3, guides = "collect")
ggsave(
    here(glue::glue("raji/output/score_stability/all_markers.pdf")),
    p_all,
    width = 15,
    height = 5
)


# Scatter 100 vs 500 permutations: one plot per cell (cell_1 = smallest … cell_5 = largest)
# cell_selection is already ordered by n_umi (ascending)
scatter_plot_list <- lapply(seq_len(5), function(i) {
    comp <- cell_selection[i]
    cell_id <- paste0("cell_", i)
    cell_size_k <- round(cell_selection_sizes[comp] / 1000, 2)
    cell_label <- paste0(cell_id, " (", cell_size_k, "k UMI)")

    scatter_data <- res_all |>
        filter(iterations %in% c("100_permutations", "500_permutations")) |>
        filter(component == comp) |>
        filter(rep == 1) |>
        pivot_wider(names_from = iterations, values_from = log2_ratio)

    cor_val <- cor(
        scatter_data$`100_permutations`,
        scatter_data$`500_permutations`,
        use = "pairwise.complete.obs"
    )

    x_min <- min(scatter_data$`100_permutations`, na.rm = TRUE)
    y_max <- max(scatter_data$`500_permutations`, na.rm = TRUE)
    library(ggrastr)
    ggplot(scatter_data, aes(x = `100_permutations`, y = `500_permutations`)) +
        geom_point_rast(size = 0.2, raster.dpi = 300) +
        geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed",
            color = "gray"
        ) +
        theme_bw(base_size = 12) +
        theme(
            aspect.ratio = 1,
            plot.title = element_text(size = 13, face = "bold"),
            axis.title = element_text(size = 11),
            axis.text = element_text(size = 10)
        ) +
        labs(
            title = cell_label,
            x = "100 permutations",
            y = "500 permutations"
        ) +
        annotate(
            "text",
            x = x_min,
            y = y_max,
            hjust = 0,
            vjust = 1,
            label = paste0("r = ", round(cor_val, 3)),
            size = 4,
            fontface = "bold"
        )
})

# Combine in one row (cell_1 … cell_5 left to right = smallest to largest)
p_scatter <- wrap_plots(scatter_plot_list, nrow = 1) +
    plot_annotation(
        title = "Correlation of scores: 100 vs 500 permutations",
        theme = theme(
            plot.title = element_text(size = 14, face = "bold")
        )
    )

ggsave(
    here("raji/output/score_stability/scatterplot_100_500_by_cell.pdf"),
    p_scatter,
    width = 14,
    height = 4
)
