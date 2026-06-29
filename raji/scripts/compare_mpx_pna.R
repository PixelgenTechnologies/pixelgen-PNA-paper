library(pixelatorR)
library(here)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)
here::i_am("raji/scripts/compare_mpx_pna.R")


########################################################
# Preprocess the MPX object
########################################################

# Read MPX data
mpx_raji <- ReadMPX_Seurat(here(
    "mpx_raji_data/Sample03_Raji_control.layout.dataset.pxl"
))

pna_raji <- ReadPNA_Seurat(here(
    "raji/data_untracked/PNA062_Raji_1000cells_S06_S6.analysis.pxl"
))

common_markers <- intersect(rownames(mpx_raji), rownames(pna_raji))

mpx_raji <- mpx_raji[common_markers, ]
pna_raji <- pna_raji[common_markers, ]

# QC MPX data
# Create a molecule rank plot to visualize the distribution of molecules per cell
# This helps identify low-quality cells with insufficient molecular counts
mpx_molecule_rank_plot <- MoleculeRankPlot(mpx_raji) +
    labs(title = "MPX Raji")

# Set a quality threshold: cells with fewer than 10,000 molecules are likely
# low-quality or damaged cells that should be excluded from downstream analysis
cutoff <- 10000
# Visualize the cutoff on the plot to verify it's appropriate for this dataset
mpx_molecule_rank_plot +
    geom_hline(yintercept = cutoff, linetype = "dashed") +
    labs(title = "MPX Raji")

# Filter cells to have at least 10000 edges
# This removes low-quality cells that may introduce noise or bias in the analysis
mpx_raji <-
    mpx_raji %>%
    subset(molecules >= cutoff)

# Create diagnostic plots to assess data quality after filtering
# Violin plots show the distribution of key QC metrics across all cells
# Log scale is used because these metrics typically span several orders of magnitude
mpx_raji[[]] %>%
    select(
        molecules,
        mean_molecules_per_a_pixel,
        reads,
        mean_reads_per_molecule
    ) %>%
    pivot_longer(
        cols = c(
            "molecules",
            "mean_molecules_per_a_pixel",
            "reads",
            "mean_reads_per_molecule"
        ),
        names_to = "metric",
        values_to = "value"
    ) %>%
    ggplot(aes(x = metric, y = value)) +
    geom_violin(draw_quantiles = 0.5, fill = "gray") +
    facet_wrap(~metric, scales = "free") +
    scale_y_log10() +
    theme_bw() +
    labs(title = "MPX Raji")

# Tau plot visualizes the spatial distribution characteristics of molecules
# This helps identify components with abnormal spatial patterns that may indicate artifacts
TauPlot(mpx_raji) + theme_bw() + theme(aspect.ratio = 1)
# Only keep the components where tau_type is normal
# Components with abnormal tau values may represent technical artifacts or
# non-biological structures that should be excluded
mpx_raji <-
    mpx_raji %>%
    subset(tau_type == "normal")

# Normalize the data using centered log-ratio (CLR) transformation
# CLR normalization is appropriate for compositional data (like marker abundances)
# as it accounts for the compositional nature and makes data more comparable across cells
mpx_raji <- Normalize(mpx_raji, method = "clr")

# Collect the isotype controls
# Isotype controls are antibodies that should not bind specifically to any target
# They serve as negative controls to estimate background/non-specific binding levels
isotype_controls <- c("mIgG1", "mIgG2a", "mIgG2b")
# Plot isotype abundance distribution to assess background signal levels
# This helps determine appropriate thresholds for filtering low-abundance markers
# that may be indistinguishable from background noise
mpx_raji %>%
    JoinLayers() %>%
    LayerData(layer = "counts") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("marker") %>%
    tidyr::pivot_longer(
        -marker,
        names_to = "component",
        values_to = "counts"
    ) %>%
    filter(marker %in% isotype_controls) %>%
    group_by(component) %>%
    summarize(iso_max = max(counts)) %>%
    ggplot(aes(iso_max)) +
    geom_histogram(binwidth = 1) +
    xlab("Maximum abundance of isotype control per component") +
    ylab("Number of components") +
    theme_bw()

########################################################
# Colocalization Analysis for MPX data
########################################################

# Create a filtered colocalization object.
# Start with all scores
# Colocalization scores measure the spatial co-occurrence of marker pairs
# which is important for understanding protein interactions and cellular organization
col_scores <- ColocalizationScores(mpx_raji)

# Collect the marker counts
# Extract raw count data for each marker in each component (cell)
# This will be used to filter colocalization scores based on marker abundance
cts <- mpx_raji %>%
    JoinLayers() %>%
    LayerData(layer = "counts") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("marker") %>%
    tidyr::pivot_longer(-marker, names_to = "component", values_to = "counts")
# Create an object with the maximum abundance of the isotype controls per component
# For each component, find the highest isotype control count as a measure of
# background signal level in that specific component
max_of_ctrls <- cts %>%
    filter(marker %in% isotype_controls) %>%
    group_by(component) %>%
    summarize(iso_max = max(counts))

# Add this info to the colocalization scores
# Enrich the colocalization scores with abundance information for both markers
# and the background (isotype control) levels to enable informed filtering
augmented_col_scores <- col_scores %>%
    left_join(cts, by = c("marker_1" = "marker", "component")) %>%
    rename(marker_1_count = counts) %>%
    left_join(cts, by = c("marker_2" = "marker", "component")) %>%
    rename(marker_2_count = counts) %>%
    left_join(max_of_ctrls) %>%
    rename(isotype_control_counts = iso_max) %>%
    rowwise() %>%
    mutate(lower_marker_count = min(marker_1_count, marker_2_count))


filtered_scores <- augmented_col_scores %>%
    rowwise() %>%
    filter(lower_marker_count > max(50, isotype_control_counts)) |>
    select(marker_1, marker_2, pearson_z) |>
    ungroup()

filtered_markers <- unique(c(
    filtered_scores$marker_1,
    filtered_scores$marker_2
))

filtered_scores <- augmented_col_scores %>%
    filter(
        (marker_1 == "CD53" & marker_2 %in% filtered_markers) |
            (marker_1 %in% filtered_markers & marker_2 == "CD53")
    ) |>
    select(marker_1, marker_2, pearson_z) |>
    ungroup() |>
    bind_rows(filtered_scores)

# Find the maximum count of rows for any marker_1/marker_2 combination
max_count <- filtered_scores %>%
    group_by(marker_1, marker_2) %>%
    summarise(count = n(), .groups = "drop") %>%
    pull(count) %>%
    max()

# Pad each marker_1/marker_2 combination with zero values until it reaches max_count
filtered_scores_padded <- filtered_scores %>%
    group_by(marker_1, marker_2) %>%
    group_modify(
        ~ {
            current_count <- nrow(.x)
            if (current_count < max_count) {
                # Add padding rows with pearson = 0
                padding_rows <- data.frame(
                    pearson_z = rep(0, max_count - current_count)
                )
                bind_rows(.x, padding_rows)
            } else {
                .x
            }
        }
    ) %>%
    ungroup()

# Replace filtered_scores with the padded version
filtered_scores <- filtered_scores_padded

mpx_hm_data <- filtered_scores %>%
    group_by(marker_1, marker_2) %>%
    summarise(
        n = n(),
        # Calculate percentage of non-zero values
        pct_of_total = 100 * sum(pearson_z != 0) / n(),
        mean = mean(pearson_z)
    )
#|>
#filter(pct_of_total >= 40)

# Create symmetric matrix by adding reverse pairs
# This ensures the heatmap shows both directions (marker_1 vs marker_2 and marker_2 vs marker_1)
mpx_hm_data_symmetric <- mpx_hm_data %>%
    bind_rows(
        mpx_hm_data %>%
            rename(marker_1 = marker_2, marker_2 = marker_1)
    ) %>%
    distinct(marker_1, marker_2, .keep_all = TRUE)

# Convert to matrix format for hierarchical clustering
# Create a matrix where rows and columns are markers, values are median pearson_z
marker_list <- unique(c(
    mpx_hm_data_symmetric$marker_1,
    mpx_hm_data_symmetric$marker_2
))
mean_df <- mpx_hm_data_symmetric %>%
    select(marker_1, marker_2, mean) %>%
    pivot_wider(
        names_from = marker_2,
        values_from = mean
    ) %>%
    as.data.frame()

# Set rownames and convert to matrix
rownames(mean_df) <- mean_df$marker_1
mean_df$marker_1 <- NULL
mean_matrix <- as.matrix(mean_df)

# Perform hierarchical clustering to determine marker order
# Use correlation-based distance to cluster markers with similar colocalization patterns
# This groups markers that have similar colocalization profiles across all other markers
cor_matrix <- cor(mean_matrix, use = "pairwise.complete.obs")
dist_matrix <- as.dist(1 - cor_matrix)
hc <- hclust(dist_matrix, method = "ward.D2")
marker_order <- hc$labels[hc$order]

# Reorder the symmetric data according to hierarchical clustering
# Reverse the marker order for both axes
mpx_hm_data_symmetric <- mpx_hm_data_symmetric %>%
    mutate(
        marker_1 = factor(marker_1, levels = marker_order),
        marker_2 = factor(marker_2, levels = rev(marker_order))
    )

# Calculate limits for capping legends using shared data
# Use percentiles to cap extreme values and improve visualization
mean_limits <- quantile(
    mpx_hm_data_symmetric$mean,
    probs = c(0.05, 0.95),
    na.rm = TRUE
)
# Calculate temporary pct_limits from MPX data only
# Will be recalculated with shared data after PNA data is ready
pct_limits <- quantile(
    mpx_hm_data_symmetric$pct_of_total,
    probs = c(0, 0.95),
    na.rm = TRUE
)

# Store MPX heatmap as a variable using shared markers
mpx_heatmap <- mpx_hm_data_symmetric %>%
    ggplot(aes(
        x = marker_1,
        y = marker_2,
        color = mean,
        # size = pct_of_total
    )) +
    geom_point(alpha = 0.8, stroke = 0, size = 6) +
    scale_color_gradient2(
        low = "#2166AC",
        mid = "#F7F7F7",
        high = "#B2182B",
        midpoint = 0,
        limits = mean_limits,
        oob = scales::squish,
        name = "Mean\nPearson",
        guide = guide_colorbar(
            barwidth = 0.8,
            barheight = 8,
            frame.colour = "black",
            ticks.colour = "black"
        )
    ) +
    # scale_size_continuous(
    #     name = "% of\nTotal Cells",
    #     range = c(2, 5),
    #     limits = pct_limits,
    #     guide = guide_legend(
    #         override.aes = list(color = "gray50", alpha = 0.8),
    #         nrow = 1
    #     )
    # ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 16 # Larger text
        ),
        axis.text.y = element_text(
            hjust = 1,
            size = 16 # Larger text
        ),
        axis.title = element_text(size = 20, face = "bold"), # Larger titles
        plot.title = element_text(
            size = 24, # Larger plot title
            face = "bold",
            hjust = 0.5,
            margin = margin(b = 15)
        ),
        legend.title = element_text(size = 16), # Larger legend titles
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = "right",
        legend.box = "vertical",
        legend.spacing = unit(0.5, "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black", linewidth = 0.8) # Add axis lines
    ) +

    labs(
        title = "MPX: Mean Pearson-Z Correlation (padded)",
        subtitle = "Marker colocalization in MPX data",
        x = "Marker 1",
        y = "Marker 2"
    ) +
    coord_fixed() +
    theme(
        axis.line.x = element_line(color = "black", linewidth = 0.8),
        axis.line.y = element_line(color = "black", linewidth = 0.8)
    )
mpx_heatmap

########################################################
# Polarization Analysis for MPX data
########################################################

polar_scores <- PolarizationScores(mpx_raji)

cts <- mpx_raji %>%
    JoinLayers() %>%
    LayerData(layer = "counts") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("marker") %>%
    tidyr::pivot_longer(-marker, names_to = "component", values_to = "counts")
# Create an object with the maximum abundance of the isotype controls per
max_of_ctrls <- cts %>%
    filter(marker %in% isotype_controls) %>%
    group_by(component) %>%
    summarize(iso_max = max(counts))

augmented_polar_scores <- polar_scores %>%
    left_join(cts, by = c("marker", "component")) %>%
    left_join(max_of_ctrls) %>%
    rename(isotype_control_counts = iso_max)

augmented_polar_scores <- augmented_polar_scores %>%
    rowwise() %>%
    filter(counts > max(50, isotype_control_counts))

# Filter polarization scores to select marker and morans_z
filtered_polar_scores <- augmented_polar_scores %>%
    select(marker, morans_z) |>
    ungroup()

# Find the maximum count of rows for any marker
max_count_polar <- filtered_polar_scores %>%
    group_by(marker) %>%
    summarise(count = n(), .groups = "drop") %>%
    pull(count) %>%
    max()

# Pad each marker combination with zero values until it reaches max_count_polar
filtered_polar_scores_padded <- filtered_polar_scores %>%
    group_by(marker) %>%
    group_modify(
        ~ {
            current_count <- nrow(.x)
            if (current_count < max_count_polar) {
                # Add padding rows with morans_z = 0
                padding_rows <- data.frame(
                    morans_z = rep(0, max_count_polar - current_count)
                )
                bind_rows(.x, padding_rows)
            } else {
                .x
            }
        }
    ) %>%
    ungroup()

# Replace filtered_polar_scores with the padded version
filtered_polar_scores <- filtered_polar_scores_padded

# Calculate mean morans_z similar to colocalization analysis
polar_hm_data <- filtered_polar_scores %>%
    group_by(marker) %>%
    summarise(
        n = n(),
        # Calculate percentage of non-zero values
        pct_of_total = 100 * sum(morans_z != 0) / n(),
        mean = mean(morans_z)
    )
#|>
#filter(pct_of_total >= 40)

########################################################
# Add polarization data to colocalization heatmap
########################################################

# Create diagonal data for polarization scores (where marker_1 == marker_2)
# Ensure the markers use the same factor levels as the colocalization data
polar_diagonal <- polar_hm_data %>%
    mutate(
        marker_1 = marker,
        marker_2 = marker,
        mean_morans_z = mean
    ) %>%
    # Add missing columns to match mpx_hm_data_symmetric structure
    mutate(
        n = n,
        pct_of_total = pct_of_total,
        mean = mean_morans_z
    ) %>%
    select(marker_1, marker_2, n, pct_of_total, mean) %>%
    # Set factor levels to match the colocalization data order
    mutate(
        marker_1 = factor(marker_1, levels = marker_order),
        marker_2 = factor(marker_2, levels = rev(marker_order))
    )

# Join polarization data with colocalization data for diagonal elements
# Use bind_rows but ensure factor levels are preserved
mpx_hm_data_symmetric_with_polar <- mpx_hm_data_symmetric %>%
    bind_rows(polar_diagonal) %>%
    # Ensure factor levels are maintained after binding
    mutate(
        marker_1 = factor(marker_1, levels = marker_order),
        marker_2 = factor(marker_2, levels = rev(marker_order))
    )

# MPX heatmap will be created after shared pct_limits are calculated
# (see below after PNA data is ready)

########################################################
# Preprocess the PNA data
########################################################

#pna_raji <- ReadPNA_Seurat(here(
#    "raji/data_untracked/PNA062_Raji_1000cells_S06_S6.analysis.pxl"
#))

pg_data <- Normalize(pna_raji, method = "clr")


DefaultAssay(pg_data) <- "PNA"
clust_scores <- ProximityScores(
    pg_data,
    add_marker_counts = TRUE,
    add_marker_proportions = TRUE
)

# define background threshold
background_threshold <- median(pg_data$isotype_fraction) +
    2 * mad(pg_data$isotype_fraction)

# summarize scores
proximity_summarized_log2 <- clust_scores |>
    FilterProximityScores(background_threshold_pct = background_threshold) |>
    SummarizeProximityScores(
        proximity_metric = "join_count_z",
        summary_stat = "mean"
    )

########################################################
# Create PNA heatmap similar to MPX heatmap
########################################################

# Get markers from MPX heatmap (using marker_order)
mpx_markers <- marker_order

# Get markers available in PNA data
pna_markers <- unique(c(
    proximity_summarized_log2$marker_1,
    proximity_summarized_log2$marker_2
))

# Find intersection of markers
common_markers_pna <- intersect(mpx_markers, pna_markers)

# Filter to only include common markers
pna_hm_data <- proximity_summarized_log2 |>
    filter(
        marker_1 %in% common_markers_pna & marker_2 %in% common_markers_pna
    ) |>
    select(marker_1, marker_2, mean_join_count_z, pct_detected) |>
    rename(mean = mean_join_count_z, pct_of_total = pct_detected) |>
    # Add n column for consistency
    mutate(n = 1) |>
    group_by(marker_1, marker_2) |>
    summarise(
        n = n(),
        pct_of_total = first(pct_of_total) * 100, # Convert to percentage
        mean = first(mean),
        .groups = "drop"
    )

# Create symmetric matrix by adding reverse pairs
pna_hm_data_symmetric <- pna_hm_data %>%
    bind_rows(
        pna_hm_data %>%
            rename(marker_1 = marker_2, marker_2 = marker_1)
    ) %>%
    distinct(marker_1, marker_2, .keep_all = TRUE)

# Perform hierarchical clustering for PNA data to determine marker order
# Convert to matrix format for hierarchical clustering
pna_marker_list <- unique(c(
    pna_hm_data_symmetric$marker_1,
    pna_hm_data_symmetric$marker_2
))
pna_mean_df <- pna_hm_data_symmetric %>%
    select(marker_1, marker_2, mean) %>%
    pivot_wider(
        names_from = marker_2,
        values_from = mean
    ) %>%
    as.data.frame()

# Set rownames and convert to matrix
rownames(pna_mean_df) <- pna_mean_df$marker_1
pna_mean_df$marker_1 <- NULL
pna_mean_matrix <- as.matrix(pna_mean_df)

# Perform hierarchical clustering to determine PNA-specific marker order
# Use correlation-based distance to cluster markers with similar colocalization patterns
pna_cor_matrix <- cor(pna_mean_matrix, use = "pairwise.complete.obs")
pna_dist_matrix <- as.dist(1 - pna_cor_matrix)
pna_hc <- hclust(pna_dist_matrix, method = "ward.D2")
pna_marker_order <- pna_hc$labels[pna_hc$order]

# Reorder the symmetric data according to PNA hierarchical clustering
pna_hm_data_symmetric <- pna_hm_data_symmetric |>
    mutate(
        marker_1 = factor(marker_1, levels = pna_marker_order),
        marker_2 = factor(marker_2, levels = rev(pna_marker_order))
        #marker_1 = factor(marker_1, levels = marker_order),
        #marker_2 = factor(marker_2, levels = rev(marker_order))
    )

# Calculate limits for capping legends (similar to MPX)
pna_mean_limits <- quantile(
    pna_hm_data_symmetric$mean,
    probs = c(0.05, 0.95),
    na.rm = TRUE
)
# Calculate shared pct_limits from both MPX and PNA data for consistent scaling
shared_pct_values <- c(
    mpx_hm_data_symmetric_with_polar$pct_of_total,
    pna_hm_data_symmetric$pct_of_total
)
pct_limits <- quantile(
    shared_pct_values,
    probs = c(0, 0.95),
    na.rm = TRUE
)
# Recreate MPX heatmap with shared pct_limits to ensure consistent scaling
mpx_heatmap <- mpx_hm_data_symmetric_with_polar %>%
    ggplot(aes(
        x = marker_1,
        y = marker_2,
        color = mean,
        size = pct_of_total
    )) +
    geom_point(alpha = 0.8, stroke = 0) +
    # Add border/outline to diagonal elements to indicate they show morans_z
    geom_point(
        data = mpx_hm_data_symmetric_with_polar %>%
            filter(marker_1 == marker_2),
        aes(x = marker_1, y = marker_2, size = pct_of_total),
        shape = 21,
        stroke = 1,
        color = "black",
        fill = NA,
        inherit.aes = FALSE
    ) +
    scale_color_gradient2(
        low = "#2166AC",
        mid = "#F7F7F7",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-3, 3),
        oob = scales::squish,
        name = "Mean Pearson Z",
        guide = guide_colorbar(
            barwidth = 8,
            barheight = 0.8,
            frame.colour = "black",
            ticks.colour = "black",
            title.position = "top",
            direction = "horizontal"
        )
    ) +
    scale_size_continuous(
        name = "% of Total Cells",
        range = c(2, 6),
        limits = pct_limits,
        guide = guide_legend(
            override.aes = list(color = "gray50", alpha = 0.8),
            direction = "horizontal",
            title.position = "top"
        )
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 14 # Larger text
        ),
        axis.text.y = element_text(
            hjust = 1,
            size = 14 # Larger text
        ),
        axis.title = element_text(size = 20, face = "bold"), # Larger titles
        plot.title = element_text(
            size = 24, # Larger plot title
            face = "bold",
            hjust = 0.5,
            margin = margin(b = 15)
        ),
        legend.title = element_text(size = 16), # Larger legend titles
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing = unit(0.5, "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black", linewidth = 0.8) # Add axis lines
    ) +

    labs(
        title = "MPX: Mean Pearson-Z Correlation",
        subtitle = "Marker colocalization in MPX data (diagonal shows mean Moran's Z)",
        x = "Marker 1",
        y = "Marker 2"
    ) +
    coord_fixed() +
    theme(
        axis.line.x = element_line(color = "black", linewidth = 0.8),
        axis.line.y = element_line(color = "black", linewidth = 0.8)
    )

# Create PNA heatmap with same style as MPX heatmap,
# adding an outlined border for diagonal elements (marker_1 == marker_2)
pna_heatmap <- pna_hm_data_symmetric %>%
    ggplot(aes(
        x = marker_1,
        y = marker_2,
        color = mean,
        size = pct_of_total
    )) +
    geom_point(alpha = 0.8, stroke = 0) +
    # Add border/outline to diagonal elements as in MPX heatmap
    geom_point(
        data = pna_hm_data_symmetric %>% filter(marker_1 == marker_2),
        aes(x = marker_1, y = marker_2, size = pct_of_total),
        shape = 21,
        stroke = 1,
        color = "black",
        fill = NA,
        inherit.aes = FALSE
    ) +
    scale_color_gradient2(
        low = "#2166AC",
        mid = "#F7F7F7",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-3, 3),
        oob = scales::squish,
        name = "Mean Join Count Z",
        guide = guide_colorbar(
            barwidth = 8,
            barheight = 0.8,
            frame.colour = "black",
            ticks.colour = "black",
            title.position = "top",
            direction = "horizontal"
        )
    ) +
    scale_size_continuous(
        name = "% of\nTotal Cells",
        range = c(2, 6),
        limits = pct_limits,
        guide = "none" # Hide size legend for PNA plot (only show on MPX)
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 14
        ),
        axis.text.y = element_text(
            hjust = 1,
            size = 14
        ),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(
            size = 24,
            face = "bold",
            hjust = 0.5,
            margin = margin(b = 15)
        ),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing = unit(0.5, "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black", linewidth = 0.8)
    ) +
    labs(
        title = "PNA: Mean Join Count Z",
        subtitle = "Marker colocalization in PNA data (diagonal outlined)",
        x = "Marker 1",
        y = "Marker 2"
    ) +
    coord_fixed() +
    theme(
        axis.line.x = element_line(color = "black", linewidth = 0.8),
        axis.line.y = element_line(color = "black", linewidth = 0.8)
    )


mpx_pna_heatmap <- mpx_heatmap + pna_heatmap

ggsave(
    here::here("output/mpx_vs_pna/mpx_vs_pna_heatmap.pdf"),
    mpx_pna_heatmap,
    width = 20,
    height = 10
)


#### Networks
markers_to_plot <- c("CD54", "CD19", "CD37", "CD82", "CD29", "CD49D")
markers_to_plot <- c("B2M", "HLA-ABC", "HLA-DR")

df_net <- ColocalizationScores(mpx_raji) %>%
    filter(marker_1 %in% markers_to_plot & marker_2 %in% markers_to_plot) %>%
    select(marker_1, marker_2, pearson_z, component) %>%
    group_by(marker_1, marker_2) %>%
    summarize(average_coloc = mean(pearson_z)) %>%
    mutate(from = marker_1, to = marker_2)

df_net_pna <- ProximityScores(pna_raji) %>%
    filter(marker_1 %in% markers_to_plot & marker_2 %in% markers_to_plot) %>%
    select(marker_1, marker_2, join_count_z) %>%
    group_by(marker_1, marker_2) %>%
    summarize(average_coloc = mean(join_count_z)) %>%
    mutate(from = marker_1, to = marker_2) %>%
    ungroup()

df_net |>
    mutate(method = "MPX") |>
    bind_rows(
        df_net_pna %>%
            mutate(method = "PNA")
    ) |>
    mutate(average_coloc = ifelse(average_coloc > 1, 1, average_coloc)) |>
    mutate(average_coloc = ifelse(average_coloc < -1, -1, average_coloc)) |>
    tidygraph::as_tbl_graph() %>%
    ggraph(layout = "kk") +
    geom_edge_link(aes(color = average_coloc), width = 2) +
    geom_node_label(
        aes(label = name),
        size = 3
    ) +
    scale_edge_color_gradientn(
        colours = c("#728BB1", "white", "#CD6F8D"),
        limits = c(-1, 1)
    ) +
    facet_edges(~method, nrow = 1) +
    theme(
        aspect.ratio = 1,
        panel.background = element_rect(fill = "white"),
        panel.spacing = unit(1, "lines")
    )


#####
# Recalculate proximity PNA wit shared markers only
########################################################
# Expand cell graphs and recompute scores in a single batched pipeline
library(tidygraph)

pna_data_subset <- pg_data


# Process in batches of cells drawn directly from pg_data
batch_size <- 10
cell_ids <- colnames(pna_data_subset)
n_cells <- length(cell_ids)
n_batches <- ceiling(n_cells / batch_size)

# Ensure the last batch has more than 1 cell to avoid validation errors
# If the last batch would only have 1 cell, merge it with the previous batch
last_batch_size <- n_cells - (n_batches - 1) * batch_size
if (last_batch_size < 2) {
    n_batches <- n_batches - 1
}

# Store per-cell score data frames (one table per cell, 991 total)
# Each element will be a data frame with proximity scores for that cell
recalculated_prox_scores_list <- vector("list", n_cells)

for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_cells)
    batch_indices <- start_idx:end_idx
    batch_cells <- cell_ids[batch_indices]

    cat(sprintf(
        "Processing batch %d/%d (cells %d-%d)...\n",
        i,
        n_batches,
        start_idx,
        end_idx
    ))

    # Restrict pg_data to this batch of cells, then build cell graphs
    pg_data_batch <- subset(pna_data_subset, cells = batch_cells)
    # Load cell graphs - this modifies pg_data_batch in place
    pg_data_batch <- LoadCellGraphs(pg_data_batch)
    # Extract cell graphs directly from the assay to avoid CellGraphs() validation error
    # (CellGraphs() fails when there's only 1 cell due to layer dimension validation)
    cg_batch <- lapply(batch_cells, function(cid) {
        pg_data_batch@assays$PNA@cellgraphs[[cid]]
    })

    # Process each cell individually to get per-cell scores
    for (j in seq_along(cg_batch)) {
        # Calculate the global cell index (across all batches)
        cell_idx <- start_idx + j - 1
        # Get the cell ID for this cell
        cell_id <- batch_cells[j]
        cg <- cg_batch[[j]]

        # Compute proximity scores for this single cell (both expanded and original)
        # ComputeProximityScores expects a named list, so name the list element with the cell_id
        # Add cell_id column to each result before storing
        recalculated_scores <- pixelatorRinternal::ComputeProximityScores(cg)
        recalculated_scores$cell_id <- cell_id
        recalculated_prox_scores_list[[cell_idx]] <- recalculated_scores
    }
}
cat("All batches processed.\n")

# Combine all per-cell scores into single result tables
# Each list contains 991 tables (one per cell), which we combine into single data frames
exp_prox_scores <- bind_rows(exp_prox_scores_list)
prox_scores_from_original <- bind_rows(prox_scores_from_original_list)


test <- exp_prox_scores |>
    filter(
        marker_1 %in% common_markers_pna & marker_2 %in% common_markers_pna
    ) |>
    select(marker_1, marker_2, metric = join_count_z) |>
    rowwise() |>
    mutate(
        # Order markers alphabetically so marker_1 is always first
        marker_1_ordered = min(marker_1, marker_2),
        marker_2_ordered = max(marker_1, marker_2)
    ) |>
    ungroup() |>
    select(marker_1 = marker_1_ordered, marker_2 = marker_2_ordered, metric)

# Find the maximum count of rows for any marker
max_count_polar_test <- test %>%
    group_by(marker_1, marker_2) %>%
    summarise(count = n(), .groups = "drop") %>%
    pull(count) %>%
    max()

# Pad each marker combination with zero values until it reaches max_count_polar
test_padded <- test %>%
    group_by(marker_1, marker_2) %>%
    group_modify(
        ~ {
            current_count <- nrow(.x)
            if (current_count < max_count_polar_test) {
                # Add padding rows with morans_z = 0
                padding_rows <- data.frame(
                    metric = rep(0, max_count_polar_test - current_count)
                )
                bind_rows(.x, padding_rows)
            } else {
                .x
            }
        }
    ) %>%
    ungroup()


test_padded_hm_data <- test_padded %>%
    group_by(marker_1, marker_2) %>%
    summarise(mean = mean(metric), .groups = "drop")


# Create symmetric matrix by adding reverse pairs
pna_hm_data_symmetric_test <- test_padded_hm_data %>%
    bind_rows(
        test_padded_hm_data %>%
            rename(marker_1 = marker_2, marker_2 = marker_1)
    ) %>%
    distinct(marker_1, marker_2, .keep_all = TRUE) |>
    # Ensure factor levels match MPX heatmap order
    mutate(
        marker_1 = factor(marker_1, levels = marker_order),
        marker_2 = factor(marker_2, levels = rev(marker_order))
    )

# Calculate limits for capping legends (similar to MPX)
pna_mean_limits <- quantile(
    pna_hm_data_symmetric_test$mean,
    probs = c(0.05, 0.95),
    na.rm = TRUE
)


# Create PNA heatmap with same style as MPX heatmap,
# adding an outlined border for diagonal elements (marker_1 == marker_2)
pna_heatmap_test <- pna_hm_data_symmetric_test %>%
    ggplot(aes(
        x = marker_1,
        y = marker_2,
        color = mean
    )) +
    geom_point(alpha = 0.8, stroke = 0, size = 6) +
    # Add border/outline to diagonal elements as in MPX heatmap
    geom_point(
        data = pna_hm_data_symmetric_test %>% filter(marker_1 == marker_2),
        aes(x = marker_1, y = marker_2),
        shape = 21,
        stroke = 1,
        color = "black",
        fill = NA,
        inherit.aes = FALSE
    ) +
    scale_color_gradient2(
        low = "#2166AC",
        mid = "#F7F7F7",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-3, 3),
        oob = scales::squish,
        name = "Mean Join Count Z",
        guide = guide_colorbar(
            barwidth = 8,
            barheight = 0.8,
            frame.colour = "black",
            ticks.colour = "black",
            title.position = "top",
            direction = "horizontal"
        )
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 14
        ),
        axis.text.y = element_text(
            hjust = 1,
            size = 14
        ),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(
            size = 24,
            face = "bold",
            hjust = 0.5,
            margin = margin(b = 15)
        ),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing = unit(0.5, "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black", linewidth = 0.8)
    ) +
    labs(
        title = "PNA: Mean Join Count Z",
        subtitle = "Marker colocalization in PNA data (diagonal outlined)",
        x = "Marker 1",
        y = "Marker 2"
    ) +
    coord_fixed() +
    theme(
        axis.line.x = element_line(color = "black", linewidth = 0.8),
        axis.line.y = element_line(color = "black", linewidth = 0.8)
    )
pna_heatmap_test


########################################################
# Compare connectivity
########################################################

# Sample 5 cells from MPX and PNA; build layouts and mid-z slices for both
set.seed(145)
mpx_sample_cells <- sample(colnames(mpx_raji), size = 5)

mpx_raji <- RemoveCellGraphs(mpx_raji)
mpx_raji <- mpx_raji %>%
    LoadCellGraphs(
        cells = mpx_sample_cells,
        load_layouts = TRUE,
        force = TRUE
    ) %>%
    ComputeLayout(layout_method = "wpmds", dim = 3)

Plot2DGraph(mpx_raji, cells = mpx_sample_cells, layout_method = "wpmds_3d")

mpx_cg_list <- CellGraphs(mpx_raji)[mpx_sample_cells]

mpx_slices <- lapply(names(mpx_cg_list), function(name) {
    CellGraphData(mpx_cg_list[[name]], "layout")$wpmds_3d %>%
        normalize_layout_coordinates() %>%
        mutate(component = name, assay = "MPX") %>%
        filter(between(z, quantile(z, 0.35), quantile(z, 0.65)))
}) %>%
    bind_rows()

set.seed(245)
pna_sample_cells <- sample(colnames(pg_data), size = 5)

pg_data <- RemoveCellGraphs(pg_data)
pg_data <- pg_data %>%
    LoadCellGraphs(
        cells = pna_sample_cells,
        load_layouts = TRUE,
        force = TRUE
    ) %>%
    ComputeLayout(layout_method = "wpmds", dim = 3)

pg_cg_list <- CellGraphs(pg_data)[pna_sample_cells]

pna_slices <- lapply(names(pg_cg_list), function(name) {
    CellGraphData(pg_cg_list[[name]], "layout")$wpmds_3d %>%
        normalize_layout_coordinates() %>%
        mutate(component = name, assay = "PNA") %>%
        filter(between(z, quantile(z, 0.35), quantile(z, 0.65)))
}) %>%
    bind_rows()

# Combine and add panel index so MPX vs PNA line up (cell 1, 2, ... 5)
xyz_combined <- bind_rows(mpx_slices, pna_slices) %>%
    group_by(assay) %>%
    mutate(
        cell_idx = as.integer(factor(component, levels = unique(component)))
    ) %>%
    ungroup()

# Single plot: MPX and PNA side by side, shared x/y scale (facet_grid scales = "fixed")
ggplot(xyz_combined, aes(x, y)) +
    geom_density_2d_filled() +
    geom_point(size = 0.2, alpha = 0.1) +
    coord_fixed() +
    theme_void() +
    facet_grid(assay ~ cell_idx) +
    scale_fill_manual(values = viridis::magma(n = 13, direction = -1))


library(tidygraph)
library(ggraph)
library(igraph)
# Compute coreness and plot its distribution

sample_cells <- sample(colnames(mpx_raji), size = 5)

mpx_raji <- RemoveCellGraphs(mpx_raji)
mpx_raji <- mpx_raji %>%
    LoadCellGraphs()
mpx_cg_list <- CellGraphs(mpx_raji)

mpx_coreness <- lapply(mpx_cg_list, function(cg) {
    cg@cellgraph %>%
        mutate(coreness = node_coreness()) %>%
        activate(nodes) %>%
        as_tibble() |>
        summarise(mean_coreness = mean(coreness), .groups = "drop")
}) |>
    bind_rows()


# Single ggplot: density of mean coreness (MPX) vs average k-core (PNA)
bind_rows(
    tibble(method = "MPX", value = mpx_coreness$mean_coreness),
    tibble(method = "PNA", value = pg_data$average_k_core)
) %>%
    ggplot(aes(x = value, color = method, fill = method)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    labs(
        x = "Coreness",
        y = "Density",
        title = "Coreness distribution: MPX vs PNA"
    ) +
    scale_color_manual(values = c(MPX = "#E41A1C", PNA = "#377EB8")) +
    scale_fill_manual(values = c(MPX = "#E41A1C", PNA = "#377EB8")) +
    theme_minimal() +
    theme(legend.title = element_blank())


########################################################
# Compare MPX and PNA heatmaps
########################################################

# Correlate the heatmapdata
mpx_heatmap_data <- mpx_heatmap$data
pna_heatmap_data <- pna_heatmap$data

# Correlate the heatmap data
mpx_heatmap_data <- mpx_heatmap_data %>%
    select(marker_1, marker_2, mean)
pna_heatmap_data <- pna_heatmap_data %>%
    select(marker_1, marker_2, mean)

# Remove symmetric duplicates by keeping canonical pairs (marker_1 <= marker_2)
mpx_heatmap_data <- mpx_heatmap_data %>%
    filter(!(marker_1 == marker_2)) %>%
    rowwise() %>%
    mutate(
        pair_1 = min(as.character(marker_1), as.character(marker_2)),
        pair_2 = max(as.character(marker_1), as.character(marker_2))
    ) %>%
    ungroup() %>%
    distinct(pair_1, pair_2, .keep_all = TRUE) %>%
    select(marker_1, marker_2, mean)

pna_heatmap_data <- pna_heatmap_data %>%
    filter(!(marker_1 == marker_2)) %>%
    rowwise() %>%
    mutate(
        pair_1 = min(as.character(marker_1), as.character(marker_2)),
        pair_2 = max(as.character(marker_1), as.character(marker_2))
    ) %>%
    ungroup() %>%
    distinct(pair_1, pair_2, .keep_all = TRUE) %>%
    select(marker_1, marker_2, mean)

library(ggrepel)

corr_coloc_data <- mpx_heatmap_data |>
    left_join(pna_heatmap_data, by = c("marker_1", "marker_2")) |>
    mutate(
        mpx_data = mean.x,
        pna_data = mean.y,
        same_sign = sign(mpx_data) == sign(pna_data) &
            !is.na(mpx_data) &
            !is.na(pna_data),
        outlier = abs(mpx_data - median(mpx_data, na.rm = TRUE)) >
            2 * sd(mpx_data, na.rm = TRUE) |
            abs(pna_data - median(pna_data, na.rm = TRUE)) >
                2 * sd(pna_data, na.rm = TRUE)
    ) |>
    ggplot(aes(x = mpx_data, y = pna_data)) +
    geom_hline(
        yintercept = 0,
        color = "gray70",
        linetype = "solid",
        linewidth = 0.5
    ) +
    geom_vline(
        xintercept = 0,
        color = "gray70",
        linetype = "solid",
        linewidth = 0.5
    ) +
    geom_point(
        alpha = 0.6,
        size = 2.5,
        color = "#2166AC",
        shape = 19
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        color = "#B2182B",
        linetype = "dashed",
        linewidth = 1
    ) +
    geom_smooth(
        method = "lm",
        se = TRUE,
        color = "#2166AC",
        fill = "#2166AC",
        alpha = 0.2,
        linewidth = 1,
        inherit.aes = FALSE,
        aes(x = mpx_data, y = pna_data)
    ) +
    xlim(-5, 5) +
    ylim(-5, 5) +
    # Label outliers with marker_1 and marker_2
    ggrepel::geom_text_repel(
        data = function(d) {
            d[d$outlier & !is.na(d$mpx_data) & !is.na(d$pna_data), ]
        },
        aes(label = paste(marker_1, marker_2, sep = " - ")),
        color = "black",
        size = 3,
        max.overlaps = 15
    ) +
    scale_color_manual(
        values = c("FALSE" = "grey60", "TRUE" = "darkgreen"),
        name = "Same Sign",
        labels = c("FALSE" = "Different", "TRUE" = "Same")
    ) +
    theme_bw() +
    theme(
        aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(
            size = 18,
            face = "bold",
            hjust = 0.5,
            margin = margin(b = 15)
        ),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed() +
    labs(
        title = "Correlation between MPX and PNA Colocalization",
        x = "Mean Pearson Z (MPX)",
        y = "Mean Join Count Z (PNA)"
    )

#### Comparison of MPX and PNA
# filtered datasets + molecules vs n_umi
molecule_comparison_data <- mpx_raji@meta.data |>
    mutate(source = "MPX") |>
    select(source, molecules, reads) |>
    bind_rows(
        pg_data@meta.data |>
            mutate(source = "PNA") |>
            select(source, molecules = n_umi, reads = reads_in_component)
    ) |>
    mutate(molecules_per_read = molecules / reads)

# Calculate medians for annotation
median_values <- molecule_comparison_data |>
    group_by(source) |>
    summarise(
        median = median(molecules_per_read, na.rm = TRUE),
        .groups = "drop"
    )

p_mol <- molecule_comparison_data |>
    ggplot(aes(x = source, y = molecules_per_read, fill = source)) +
    geom_violin() +
    stat_summary(
        fun = median,
        geom = "point",
        size = 3,
        color = "black",
        shape = 21,
        fill = "white",
        stroke = 1.5
    ) +
    geom_text(
        data = median_values,
        aes(x = source, y = median, label = round(median, 2)),
        vjust = -1.5,
        size = 5,
        fontface = "bold",
        inherit.aes = FALSE
    ) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(
            size = 16,
            face = "bold",
            hjust = 0.5,
            vjust = -1.5,
            margin = margin(b = 15)
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
    ) +
    theme(aspect.ratio = 1) +
    labs(
        title = "",
        x = "",
        y = "Number of molecules per read"
    ) +
    ylim(0.2, 0.6)


# filtered datasets + molecules vs n_umi
antibody_comparison_data <- FetchData(mpx_raji, vars = rownames(mpx_raji)) |>
    tibble::rownames_to_column("component") |>
    pivot_longer(
        cols = -component,
        names_to = "marker",
        values_to = "expression"
    ) |>
    filter(expression > 1) |>
    group_by(component) |>
    summarise(n_markers = n(), .groups = "drop") |>
    mutate(source = "MPX") |>
    bind_rows(
        FetchData(pg_data, vars = rownames(pg_data)) |>
            tibble::rownames_to_column("component") |>
            pivot_longer(
                cols = -component,
                names_to = "marker",
                values_to = "expression"
            ) |>
            filter(expression > 1) |>
            group_by(component) |>
            summarise(n_markers = n(), .groups = "drop") |>
            mutate(source = "PNA")
    )

# Calculate medians for annotation
median_values <- antibody_comparison_data |>
    group_by(source) |>
    summarise(median = median(n_markers, na.rm = TRUE), .groups = "drop")

p_markers <- antibody_comparison_data |>
    ggplot(aes(x = source, y = n_markers, fill = source)) +
    geom_violin() +
    stat_summary(
        fun = median,
        geom = "point",
        size = 3,
        color = "black",
        shape = 21,
        fill = "white",
        stroke = 1.5
    ) +
    geom_text(
        data = median_values,
        aes(x = source, y = median, label = round(median, 0)),
        vjust = -1.5,
        size = 5,
        fontface = "bold",
        inherit.aes = FALSE
    ) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(
            size = 16,
            face = "bold",
            hjust = 0.5,
            vjust = -1.5,
            margin = margin(b = 15)
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
    ) +
    theme(aspect.ratio = 1) +
    labs(
        title = "",
        x = "",
        y = "Number of Markers (CLR > 1)"
    )

p_mol + p_markers

ggsave(
    here::here("output/mpx_vs_pna/molecule_vs_antibody_comparison.pdf"),
    p_mol / p_markers,
    width = 20,
    height = 10
)

# Compare Abundance of markers in MPX and PNA
abundance_comparison_data <- mpx_raji@assays$mpxCells$data |>
    as.data.frame() |>
    tibble::rownames_to_column("marker") |>
    pivot_longer(
        cols = -marker,
        names_to = "component",
        values_to = "abundance"
    ) |>
    group_by(marker) |>
    summarise(
        mpx_abundance = median(abundance, na.rm = TRUE),
        .groups = "drop"
    ) |>
    full_join(
        pg_data@assays$PNA$data |>
            as.data.frame() |>
            tibble::rownames_to_column("marker") |>
            pivot_longer(
                cols = -marker,
                names_to = "component",
                values_to = "abundance"
            ) |>
            group_by(marker) |>
            summarise(
                pna_abundance = median(abundance, na.rm = TRUE),
                .groups = "drop"
            ),
        by = "marker"
    )

# Identify the 6 most abundant markers (using combined abundance score)
top_markers <- abundance_comparison_data |>
    filter(!is.na(mpx_abundance) & !is.na(pna_abundance)) |>
    mutate(
        combined_abundance = (mpx_abundance + pna_abundance) / 2
    ) |>
    slice_max(combined_abundance, n = 6) |>
    pull(marker)

p_abundance_plot <- abundance_comparison_data |>
    ggplot(aes(x = mpx_abundance, y = pna_abundance)) +
    geom_hline(
        yintercept = 0,
        color = "gray70",
        linetype = "solid",
        linewidth = 0.5
    ) +
    geom_vline(
        xintercept = 0,
        color = "gray70",
        linetype = "solid",
        linewidth = 0.5
    ) +
    geom_point(
        alpha = 0.6,
        size = 2.5,
        color = "#2166AC",
        shape = 19
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        color = "#B2182B",
        linetype = "dashed",
        linewidth = 1
    ) +
    geom_smooth(
        method = "lm",
        se = TRUE,
        color = "#2166AC",
        fill = "#2166AC",
        alpha = 0.2,
        linewidth = 1
    ) +
    ggrepel::geom_text_repel(
        data = abundance_comparison_data |>
            filter(
                marker %in%
                    top_markers &
                    !is.na(mpx_abundance) &
                    !is.na(pna_abundance)
            ),
        aes(label = marker),
        color = "black",
        size = 3.5,
        fontface = "bold",
        max.overlaps = 10,
        box.padding = 0.5,
        point.padding = 0.3
    ) +
    theme_bw() +
    theme(
        aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(
            size = 18,
            face = "bold",
            hjust = 0.5,
            margin = margin(b = 15)
        ),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed() +
    labs(
        title = "Correlation between MPX and PNA Abundance [shared markers only]",
        x = "MPX Abundance (median CLR)",
        y = "PNA Abundance (median CLR)"
    ) +
    xlim(-3, 7) +
    ylim(-3, 7)

panel_plot <- mpx_pna_heatmap /
    (p_mol + p_markers) |
    (corr_coloc_data / p_abundance_plot)
ggsave(
    here::here("output/mpx_vs_pna/panel_plot.pdf"),
    panel_plot,
    width = 30,
    height = 20
)

### COMPARE SCORE DISTRIBUTIONS
mpx_colocalization_scores <- ColocalizationScores(mpx_raji)
mpx_polarization_scores <- PolarizationScores(mpx_raji)
pna_colocalization_scores <- ProximityScores(
    pg_data,
    add_marker_counts = TRUE,
    add_marker_proportions = TRUE
)

# Correlate the heatmap data
mpx_heatmap_data <- mpx_heatmap$data
pna_heatmap_data <- pna_heatmap$data

shared_markers <- as.vector(unique(c(
    mpx_heatmap_data$marker_1,
    mpx_heatmap_data$marker_2,
    pna_heatmap_data$marker_1,
    pna_heatmap_data$marker_2
)))

mpx_colocalization_scores_filtered_pearson_z <- mpx_colocalization_scores %>%
    filter(marker_1 %in% shared_markers & marker_2 %in% shared_markers) |>
    select(marker_1, marker_2, metric = pearson_z) |>
    mutate(source = "MPX - pearson Z")

mpx_polarization_scores_filtered_morans_z <- mpx_polarization_scores %>%
    filter(marker %in% shared_markers) |>
    select(marker, metric = morans_z) |>
    mutate(marker_1 = marker, marker_2 = marker, source = "MPX - morans Z") |>
    select(-marker)

mpx_scores_filtered <- bind_rows(
    mpx_colocalization_scores_filtered_pearson_z,
    mpx_polarization_scores_filtered_morans_z
)

pna_colocalization_scores_filtered_log2_ratio <- pna_colocalization_scores %>%
    filter(marker_1 %in% shared_markers & marker_2 %in% shared_markers) |>
    select(marker_1, marker_2, metric = log2_ratio) |>
    mutate(source = "PNA - log2 ratio")

pna_colocalization_scores_filtered_join_count_z <- pna_colocalization_scores %>%
    filter(marker_1 %in% shared_markers & marker_2 %in% shared_markers) |>
    select(marker_1, marker_2, metric = join_count_z) |>
    mutate(source = "PNA - join count Z")

raw_scores_all <- bind_rows(
    mpx_scores_filtered,
    pna_colocalization_scores_filtered_log2_ratio,
    pna_colocalization_scores_filtered_join_count_z
)

raw_scores_all %>%
    mutate(pair = paste(marker_1, marker_2, sep = "_")) |>
    filter(source != "PNA - join count Z") |>
    # filter(pair %in% sample(unique(pair), size = 10, replace = FALSE)) |>
    filter(
        pair %in%
            c(
                "CD11a_CD18",
                "CD29_CD49D",
                "CD20_CD82",
                "CD54_CD54RB",
                "CD48_CD55",
                "CD19_CD20",
                "CD19_HLA-ABC",
                "CD71_CD78",
                "CD82_CD82",
                "CD20_CD20"
            )
    ) |>
    group_by(pair) |>
    ggplot(aes(x = pair, y = metric, fill = source)) +
    geom_boxplot(
        position = position_dodge(width = 0.8),
        width = 0.6,
        alpha = 0.7,
        outlier.alpha = 0.5
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    ) +
    ylim(-10, 10) +
    labs(
        title = "Raw Scores",
        x = "Marker Pair",
        y = "Metric Value",
        fill = "Method"
    )

# Expand cell graphs and recompute scores in a single batched pipeline
library(tidygraph)

# Process in batches of cells drawn directly from pg_data
batch_size <- 10
cell_ids <- colnames(pg_data)
n_cells <- length(cell_ids)
n_batches <- ceiling(n_cells / batch_size)

# Ensure the last batch has more than 1 cell to avoid validation errors
# If the last batch would only have 1 cell, merge it with the previous batch
last_batch_size <- n_cells - (n_batches - 1) * batch_size
if (last_batch_size < 2) {
    n_batches <- n_batches - 1
}

# Store per-cell score data frames (one table per cell, 991 total)
# Each element will be a data frame with proximity scores for that cell
exp_prox_scores_list <- vector("list", n_cells)
prox_scores_from_original_list <- vector("list", n_cells)

for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_cells)
    batch_indices <- start_idx:end_idx
    batch_cells <- cell_ids[batch_indices]

    cat(sprintf(
        "Processing batch %d/%d (cells %d-%d)...\n",
        i,
        n_batches,
        start_idx,
        end_idx
    ))

    # Restrict pg_data to this batch of cells, then build cell graphs
    pg_data_batch <- subset(pg_data, cells = batch_cells)
    # Load cell graphs - this modifies pg_data_batch in place
    pg_data_batch <- LoadCellGraphs(pg_data_batch)
    # Extract cell graphs directly from the assay to avoid CellGraphs() validation error
    # (CellGraphs() fails when there's only 1 cell due to layer dimension validation)
    cg_batch <- lapply(batch_cells, function(cid) {
        pg_data_batch@assays$PNA@cellgraphs[[cid]]
    })

    # Process each cell individually to get per-cell scores
    for (j in seq_along(cg_batch)) {
        # Calculate the global cell index (across all batches)
        cell_idx <- start_idx + j - 1
        # Get the cell ID for this cell
        cell_id <- batch_cells[j]
        cg <- cg_batch[[j]]

        # Expand this cell's graph
        # Convert cell graph to adjacency matrix
        A <- igraph::as_adjacency_matrix(cg@cellgraph, sparse = TRUE)
        # Expand adjacency matrix
        A_exp <- expand_adjacency_matrix(A, k = 3)
        # Convert expanded adjacency matrix back to tbl_graph
        tbl_graph_exp <- igraph::graph_from_adjacency_matrix(
            A_exp,
            mode = "undirected"
        ) %>%
            tidygraph::as_tbl_graph(directed = FALSE)
        # Add node type to tbl_graph
        tbl_graph_exp <- tbl_graph_exp %>%
            activate(nodes) %>%
            mutate(node_type = gsub(".*-", "", name))
        # Add tbl_graph to cell graph
        cg_data_exp <- cg
        cg_data_exp@cellgraph <- tbl_graph_exp

        # Compute proximity scores for this single cell (both expanded and original)
        # ComputeProximityScores expects a named list, so name the list element with the cell_id
        # Add cell_id column to each result before storing
        cg_list_exp <- list(cg_data_exp)
        names(cg_list_exp) <- cell_id
        exp_scores <- pixelatorRinternal::ComputeProximityScores(cg_list_exp)
        exp_scores$cell_id <- cell_id
        exp_prox_scores_list[[cell_idx]] <- exp_scores

        cg_list_orig <- list(cg)
        names(cg_list_orig) <- cell_id
        orig_scores <- pixelatorRinternal::ComputeProximityScores(cg_list_orig)
        orig_scores$cell_id <- cell_id
        prox_scores_from_original_list[[cell_idx]] <- orig_scores
    }
}
cat("All batches processed.\n")

# Combine all per-cell scores into single result tables
# Each list contains 991 tables (one per cell), which we combine into single data frames
exp_prox_scores <- bind_rows(exp_prox_scores_list)
prox_scores_from_original <- bind_rows(prox_scores_from_original_list)


test <- exp_prox_scores |>
    filter(
        marker_1 %in% common_markers_pna & marker_2 %in% common_markers_pna
    ) |>
    select(marker_1, marker_2, metric = join_count_z) |>
    rowwise() |>
    mutate(
        # Order markers alphabetically so marker_1 is always first
        marker_1_ordered = min(marker_1, marker_2),
        marker_2_ordered = max(marker_1, marker_2)
    ) |>
    ungroup() |>
    select(marker_1 = marker_1_ordered, marker_2 = marker_2_ordered, metric)

# Find the maximum count of rows for any marker
max_count_polar_test <- test %>%
    group_by(marker_1, marker_2) %>%
    summarise(count = n(), .groups = "drop") %>%
    pull(count) %>%
    max()

# Pad each marker combination with zero values until it reaches max_count_polar
test_padded <- test %>%
    group_by(marker_1, marker_2) %>%
    group_modify(
        ~ {
            current_count <- nrow(.x)
            if (current_count < max_count_polar_test) {
                # Add padding rows with morans_z = 0
                padding_rows <- data.frame(
                    metric = rep(0, max_count_polar_test - current_count)
                )
                bind_rows(.x, padding_rows)
            } else {
                .x
            }
        }
    ) %>%
    ungroup()


test_padded_hm_data <- test_padded %>%
    group_by(marker_1, marker_2) %>%
    summarise(mean = mean(metric), .groups = "drop")


# Create symmetric matrix by adding reverse pairs
pna_hm_data_symmetric_test <- test_padded_hm_data %>%
    bind_rows(
        test_padded_hm_data %>%
            rename(marker_1 = marker_2, marker_2 = marker_1)
    ) %>%
    distinct(marker_1, marker_2, .keep_all = TRUE) |>
    # Ensure factor levels match MPX heatmap order
    mutate(
        marker_1 = factor(marker_1, levels = marker_order),
        marker_2 = factor(marker_2, levels = rev(marker_order))
    )

# Calculate limits for capping legends (similar to MPX)
pna_mean_limits <- quantile(
    pna_hm_data_symmetric_test$mean,
    probs = c(0.05, 0.95),
    na.rm = TRUE
)


# Create PNA heatmap with same style as MPX heatmap,
# adding an outlined border for diagonal elements (marker_1 == marker_2)
pna_heatmap_test <- pna_hm_data_symmetric_test %>%
    ggplot(aes(
        x = marker_1,
        y = marker_2,
        color = mean
    )) +
    geom_point(alpha = 0.8, stroke = 0, size = 6) +
    # Add border/outline to diagonal elements as in MPX heatmap
    geom_point(
        data = pna_hm_data_symmetric_test %>% filter(marker_1 == marker_2),
        aes(x = marker_1, y = marker_2),
        shape = 21,
        stroke = 1,
        color = "black",
        fill = NA,
        inherit.aes = FALSE
    ) +
    scale_color_gradient2(
        low = "#2166AC",
        mid = "#F7F7F7",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-3, 3),
        oob = scales::squish,
        name = "Mean Join Count Z",
        guide = guide_colorbar(
            barwidth = 8,
            barheight = 0.8,
            frame.colour = "black",
            ticks.colour = "black",
            title.position = "top",
            direction = "horizontal"
        )
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 14
        ),
        axis.text.y = element_text(
            hjust = 1,
            size = 14
        ),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(
            size = 24,
            face = "bold",
            hjust = 0.5,
            margin = margin(b = 15)
        ),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing = unit(0.5, "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black", linewidth = 0.8)
    ) +
    labs(
        title = "PNA: Mean Join Count Z",
        subtitle = "Marker colocalization in PNA data (diagonal outlined)",
        x = "Marker 1",
        y = "Marker 2"
    ) +
    coord_fixed() +
    theme(
        axis.line.x = element_line(color = "black", linewidth = 0.8),
        axis.line.y = element_line(color = "black", linewidth = 0.8)
    )
pna_heatmap_test


### Compare connectivity

### Reviewer request: Linker distribution

#| echo: false
#| code-fold: true

## Use not subsetted
# Select 100 sample cells - Full datasets
set.seed(120)
sample_cells <- sample(colnames(mpx_raji), size = 100)
mpx_raji_sampled <- subset(mpx_raji, cells = sample_cells)
# Load cell graphs
mpx_raji_sampled <- LoadCellGraphs(mpx_raji_sampled)
mpx_cg_list <- CellGraphs(mpx_raji_sampled)


sample_cells <- sample(colnames(pna_raji), size = 100)
pna_raji_sampled <- subset(pna_raji, cells = sample_cells)
# Load cell graphs
pna_raji_sampled <- LoadCellGraphs(pna_raji_sampled)
pna_cg_list <- CellGraphs(pna_raji_sampled)


# Compute linker scores
calculate_degree <- function(cg) {
    cg@cellgraph |>
        activate(nodes) |>
        mutate(degree = centrality_degree()) |>
        as_tibble()
}

degree_list_mpx <- lapply(mpx_cg_list, calculate_degree)
degree_df_mpx <- bind_rows(degree_list_mpx, .id = "sample") |>
    mutate(degree_capped = pmin(degree, 10)) |>
    group_by(degree_capped) |>
    summarise(count = n()) |>
    mutate(percentage = count / sum(count) * 100)

degree_list_pna <- lapply(pna_cg_list, calculate_degree)
degree_df_pna <- bind_rows(degree_list_pna, .id = "sample") |>
    mutate(degree_capped = pmin(degree, 10)) |>
    group_by(degree_capped) |>
    summarise(count = n()) |>
    mutate(percentage = count / sum(count) * 100)

mean_degree_mpx <- mean(bind_rows(degree_list_mpx)$degree)
mean_degree_pna <- mean(bind_rows(degree_list_pna)$degree)

max_degree_mpx <- max(bind_rows(degree_list_mpx)$degree)
max_degree_pna <- max(bind_rows(degree_list_pna)$degree)

# Vertical line positions (cap at 20 so lines stay on plot)
x_mpx <- min(mean_degree_mpx, 10)
x_pna <- min(mean_degree_pna, 10)

data_to_plot <- bind_rows(
    degree_df_mpx |> mutate(method = "MPX"),
    degree_df_pna |> mutate(method = "PNA")
)

(bin_plot_mpx <- ggplot(
    data_to_plot,
    aes(x = factor(degree_capped), y = percentage, fill = method)
) +
    geom_bar(position = "dodge", stat = "identity", alpha = 0.8) +
    geom_vline(
        xintercept = as.numeric(factor(
            x_mpx,
            levels = sort(unique(data_to_plot$degree_capped))
        )),
        linetype = "dashed",
        color = "#E41A1C"
    ) +
    geom_vline(
        xintercept = as.numeric(factor(
            x_pna,
            levels = sort(unique(data_to_plot$degree_capped))
        )),
        linetype = "dashed",
        color = "#377EB8"
    ) +
    annotate(
        "text",
        x = as.numeric(factor(
            x_mpx,
            levels = sort(unique(data_to_plot$degree_capped))
        )),
        y = 2.8,
        label = sprintf("Mean MPX = %.1f", mean_degree_mpx),
        vjust = 0,
        hjust = -0.05,
        size = 3.5,
        color = "#E41A1C"
    ) +
    annotate(
        "text",
        x = as.numeric(factor(
            x_pna,
            levels = sort(unique(data_to_plot$degree_capped))
        )),
        y = 2.8,
        label = sprintf("Mean PNA = %.1f", mean_degree_pna),
        vjust = 0,
        hjust = 1.05,
        size = 3.5,
        color = "#377EB8"
    ) +
    annotate(
        "text",
        x = as.numeric(factor(
            x_mpx,
            levels = sort(unique(data_to_plot$degree_capped))
        )),
        y = 2.7,
        label = sprintf("Max MPX = %.1f", max_degree_mpx),
        vjust = 0,
        hjust = -0.05,
        size = 3.5,
        color = "#E41A1C"
    ) +
    annotate(
        "text",
        x = as.numeric(factor(
            x_pna,
            levels = sort(unique(data_to_plot$degree_capped))
        )),
        y = 2.7,
        label = sprintf("Max PNA = %.1f", max_degree_pna),
        vjust = 0,
        hjust = 1.05,
        size = 3.5,
        color = "#377EB8"
    ) +
    scale_x_discrete(
        limits = as.character(seq(0, degree_cap, by = 1)),
        breaks = as.character(seq(0, degree_cap, by = 5)),
        name = "Degree (values \u2265 10 binned at 10)"
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.15)),
        name = "Percentage of total (%)"
    ) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    labs(
        title = "Degree distribution (100 Raji cells sampled)",
        fill = "Method"
    ))

ggsave(
    here("raji/output/linker_distribution/degree_distribution_all.pdf"),
    bin_plot_mpx,
    width = 10,
    height = 10
)
