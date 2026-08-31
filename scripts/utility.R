
order_cd_markers <-
  function(markers, control_markers) {

    cd_markers <- str_detect(markers, "^CD\\d")
    control_markers <- markers %in% control_markers


    cd_order <-
      tibble(marker = markers[cd_markers]) %>%
      mutate(marker_i = str_remove(marker, "^CD") %>%
               str_extract("^\\d*") %>%
               as.numeric()) %>%
      arrange(marker_i, marker) %>%
      pull(marker)

    non_cd_order <-
      tibble(marker = markers[!cd_markers & !control_markers]) %>%
      arrange(marker) %>%
      pull(marker)

    control_order <-
      tibble(marker = markers[control_markers]) %>%
      arrange(marker) %>%
      pull(marker)

    return(c(cd_order, non_cd_order, control_order))
  }


plot_embedding <-
  function(object,
           plot_reduction = "pca",
           dims = 1:2,
           metavars = "seurat_clusters",
           pal = NULL,
           label = TRUE,
           plot_title = "",
           legend_position = "none",
           legend_cols = 2,
           ...) {



    embeddings <-
      Embeddings(object, reduction = plot_reduction) %>%
      as_tibble(rownames = "sample_component") %>%
      select(1, dims + 1)

    plot_metadata <-
      FetchData(object, vars = metavars) %>%
      as_tibble(rownames = "sample_component")

    plot_dimnames <-
      colnames(embeddings)[dims + 1]

    plot_data <-
      left_join(embeddings, plot_metadata,
                by = "sample_component") %>%
      rename(V1 = !!plot_dimnames[1], V2 = !!plot_dimnames[2],
             meta = !!metavars)

    p <-
      plot_data %>%
      ggplot(aes(x = V1, y = V2, color = meta)) +
      geom_point(size = 0.75) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = legend_position,
            legend.key.size = unit(0.5, "cm"),
            legend.spacing.x = unit(0.2, "cm"),
            legend.direction = "vertical") +
      coord_fixed() +
      labs(title = plot_title,
           x = plot_dimnames[1],
           y = plot_dimnames[2],
           color = "") +
      guides(color = guide_legend(
        ncol = legend_cols))

    if(!is.null(pal)) {
      p <- p + scale_color_manual(values = pal)
    }

    if(label) {
      p <-
        p +
        geom_text(data = . %>%
                    group_by(meta) %>%
                    summarise(V1 = median(V1), V2 = median(V2)),
                  aes(label = meta),
                  color = "black",
                  alpha = 0.75)

    }

    return(p)
  }

# Function to create smoothed visualization layout for one or more markers
create_marker_plot <- function(
    df_rotated,
    marker_name,
    k = k,
    cutoff = cutoff,
    order = FALSE
) {
    # Extract the marker data and prepare for neighbor calculation by coalescing the marker names into a single column called patch
    cell <- df_rotated %>%
        mutate(patch = coalesce(!!!syms(marker_name))) %>%
        select(x_rotated, y_rotated, z_rotated, patch)

    # Calculate neighbor patch values
    nn_df <- calculate_neighbor_patch_values(cell, k = k)
    df_rotated$avg_neighbor_patch <- nn_df$avg_neighbor_patch

    # Create the plot
    plot <- df_rotated %>%
        # Optionally order points by neighbor patch value to control overlapping points
        {
            if (order) arrange(., avg_neighbor_patch) else .
        } %>%
        # Filter out values below cutoff to reduce noise
        mutate(
            avg_neighbor_patch = ifelse(
                avg_neighbor_patch > cutoff,
                avg_neighbor_patch,
                0
            )
        ) %>%
        ggplot(aes(
            x = x_rotated,
            y = y_rotated,
            color = avg_neighbor_patch,
            fill = avg_neighbor_patch
        )) +
        geom_point(alpha = 0.4) +
        scale_color_gradientn(colors = c("grey", "#2955AE", "#1C3A76")) +
        scale_fill_gradientn(colors = c("grey", "#2955AE", "#1C3A76")) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank()
        ) +
        ggtitle(paste(marker_name, collapse = ", ")) +
        NoLegend()

    return(plot)
}

# Function to calculate average patch values of nearest neighbors
calculate_neighbor_patch_values <- function(
    cell_data,
    k = k,
    coords_cols = c("x_rotated", "y_rotated", "z_rotated"),
    patch_col = "patch"
) {
    # Create a matrix of just the coordinates for distance calculation
    coords <- cell_data %>%
        select(all_of(coords_cols)) %>%
        as.matrix()

    # Using the FNN package for efficient nearest neighbor search
    if (!requireNamespace("FNN", quietly = TRUE)) {
        install.packages("FNN")
        library(FNN)
    } else {
        library(FNN)
    }

    # Get indices of k nearest neighbors (includes the point itself as first neighbor)
    nn_indices <- get.knn(
        coords,
        k = k,
        algorithm = "kd_tree"
    )$nn.index

    # Convert to a data frame for easier analysis
    nn_df <- data.frame(
        point_id = 1:nrow(cell_data),
        nn_indices
    )

    # Helper function to get the average patch value for a set of neighbor indices
    get_avg_patch <- function(neighbors, patch_values) {
        mean(patch_values[neighbors], na.rm = TRUE)
    }

    # Apply the function to each row to calculate average patch value of neighbors
    nn_df$avg_neighbor_patch <- apply(
        nn_df[, 1:k],
        1,
        function(neighbors) {
            get_avg_patch(neighbors, cell_data[[patch_col]])
        }
    )

    return(nn_df)
}


#' Detect patches in a PNA graph
#'
#' A patch is defined as a group of nodes that are at least weakly connected to
#' each other and are enriched for a set of patch-specific markers. The algorithm
#' works as follows:
#' 1. Expand the graph adjacency matrix to consider the k-th neighbors.
#' 2. Subset the expanded adjacency matrix to include nodes that are labeled by the patch-specific markers as well
#' as the nodes connecting these patch-specific marker nodes.
#' 3. Create a new graph from the subsetted adjacency matrix and split the graph into its connected components.
#' 4. Run an iterative community detection (Leiden) to split up weakly connected patch components.
#' 5. Filter out connected patches with fewer than \code{patch_nodes_threshold} nodes.
#' 6. Label nodes in the original graph with the patch information. The largest "patch" is labeled as 1
#' and should typically correspond to the membrane graph. The rest of the patches are labeled as 2, 3, etc. and
#' correspond to patches ordered by decreasing size.
#'
#' @param cg A \code{CellGraph} object with PNA data.
#' @param patch_markers A character vector with the names of the markers that are
#' exclusively found on the patches, e.g. "CD41" for platelets.
#' @param k The number of steps to consider in the adjacency matrix. This will help
#' finding patches even if the patch markers are not directly connected to each other.
#' @param leiden_resolution The resolution parameter for the Leiden algorithm.
#' @param patch_nodes_threshold The minimum number of nodes to consider a patch.
#' @param verbose A logical indicating if messages should be printed to the console.
#'
#' @return A \code{CellGraph} object with the patches detected.
#'
patch_detection <- function(
    cg,
    patch_markers = c("CD41", "CD9", "CD62P"),
    k = 2L,
    leiden_resolution = 0.02,
    patch_nodes_threshold = 20,
    verbose = TRUE
) {
    if (verbose) {
        cli::cli_alert_info("Extracting connected patches...")
    }

    # Get adjacency matrix
    A <- cg@cellgraph %>% igraph::as_adjacency_matrix()
    A_exp <- pixelatorRinternal::expand_adjacency_matrix(A, k = k)
    patch_nodes <- (cg@counts[, patch_markers, drop = FALSE] %>%
        Matrix::rowSums()) ==
        1

    # Find connecting nodes
    A_d <- A_exp - A
    interconnected_nodes <- (A_d[patch_nodes, ] %>% Matrix::colSums()) > 1

    Ap <- A_exp[
        patch_nodes | interconnected_nodes,
        patch_nodes | interconnected_nodes
    ]

    if (verbose) {
        cli::cli_alert_info(
            "{.val {ncol(Ap)}} out of {.val {ncol(A)}} nodes are labelled as patch nodes..."
        )
    }

    g_patch_list <- igraph::graph_from_adjacency_matrix(
        Ap,
        mode = "undirected"
    ) %>%
        as_tbl_graph(directed = FALSE) %>%
        mutate(comp = igraph::components(.)$membership) %>%
        group_by(comp) %>%
        mutate(n = n()) %>%
        filter(n >= patch_nodes_threshold) %>%
        to_components()

    if (length(g_patch_list) == 0) {
        cli::cli_alert_info("No connected patches found...")
    }

    if (verbose) {
        cli::cli_alert_info(
            "Found {.val {length(g_patch_list)}} connected patches after filtering..."
        )
    }
    if (verbose) {
        cli::cli_alert_info(
            "Splitting up large patches using Leiden with resolution={.val {leiden_resolution}}..."
        )
    }

    g_patch_nodes <- lapply(seq_along(g_patch_list), function(i) {
        g_patch_list[[i]] %>%
            iterative_leiden(
                resolution = leiden_resolution,
                verbose = FALSE
            ) %N>%
            as_tibble() %>%
            mutate(patch = paste0(community, "_", i))
    }) %>%
        bind_rows() %>%
        select(-contains("community")) %>%
        group_by(patch) %>%
        mutate(n = n()) %>%
        filter(n >= patch_nodes_threshold)

    if (verbose) {
        cli::cli_alert_info(
            "Found {.val {length(unique(g_patch_nodes$patch))}} patches after splitting..."
        )
    }

    g_patch_nodes <- g_patch_nodes %>%
        group_by(patch) %>%
        count(name = "size") %>%
        arrange(-size) %>%
        ungroup() %>%
        mutate(group = row_number() + 1) %>%
        right_join(g_patch_nodes, by = "patch") %>%
        select(-size) %>%
        ungroup() %>%
        select(name, patch = group)

    # Add group labels to original graph
    g <- cg@cellgraph %N>%
        select(-matches("patch")) %>%
        left_join(g_patch_nodes, by = c("name" = "name")) %>%
        mutate(
            patch = case_when(
                is.na(patch) ~ 1,
                TRUE ~ patch
            )
        )

    if (verbose) {
        cli::cli_alert_success("Finished!")
    }

    cg@cellgraph <- g
    return(cg)
}


# Create a function to symmetrize a matrix
ultosymmetric <- function(m) {
    m <- m + t(m) - diag(diag(m))
    return(m)
}


# Create a function to generate proximity score heatmaps
#' Generate Proximity Score Heatmaps
#'
#' This function creates heatmaps visualizing the proximity scores between markers in specific cell types and conditions.
#'
#' @param prox_scores_summary A dataframe containing proximity scores between marker pairs
#' @param seurat_obj A Seurat object containing expression data
#' @param cell_type Character string specifying the cell type to analyze (default: "NK")
#' @param sample Character string specifying the sample condition to analyze (default: "pbmc_skbr3")
#' @param expr_threshold Numeric threshold for marker expression filtering (default: 25)
#' @param limit_values Numeric value to cap the heatmap color scale (default: 3)
#' @param metric Character string specifying the metric to use for the heatmap (default: "median_prox")
#' @param metric_name Character string specifying the name of the metric to use for the heatmap (default: "Median Log2 Ratio")
#' @param denoising Character string specifying the denoised data to use for the heatmap (default: "full")
#' @param return_order Logical indicating whether to return marker order instead of plot (default: FALSE)
#' @param set_order Optional vector of marker names to set a specific order (default: NULL). If NULL, the markers are ordered using hierarchical clustering. If not NULL, the markers are ordered as specified in set_order but subsetted to the markers present in the data.
#'
#' @return Either a ggplot object containing the heatmap or a vector of marker names if return_order=TRUE
create_proximity_heatmap <- function(
    prox_scores_summary,
    seurat_obj,
    metric = "median_prox",
    cell_type_column = "manual_annotation_broad",
    cell_type = "NK",
    condition_name = "pbmc_skbr3",
    denoising = "full",
    expr_threshold = 25,
    limit_values = 3,
    return_order = FALSE,
    set_order = NULL
) {
    metric_name <- switch(
        metric,
        "median_prox" = "Median Log2 Ratio",
        "mean_prox" = "Mean Log2 Ratio"
    )

    # Create matrix
    m <- prox_scores_summary %>%
        ungroup() %>%
        filter(!!sym(cell_type_column) == cell_type) %>%
        filter(condition == condition_name) %>%
        filter(type == denoising) %>%
        select(marker_1, marker_2, metric) %>%
        pivot_wider(
            names_from = marker_2,
            values_from = metric,
            values_fill = 0
        ) %>%
        column_to_rownames("marker_1") %>%
        as.matrix()
    m <- ultosymmetric(m)

    # Get median values per protein
    cts <- FetchData(
        seurat_obj,
        vars = c(cell_type_column, "condition", "type", rownames(seurat_obj)),
        layer = "counts"
    ) %>%
        group_by(!!sym(cell_type_column), condition, type) %>%
        summarise(across(all_of(rownames(seurat_obj)), median)) %>%
        pivot_longer(
            cols = -c(!!sym(cell_type_column), condition, type),
            names_to = "marker",
            values_to = "median_cts"
        )

    # Filter for expressed markers
    if (is.null(set_order)) {
        expr_mrk <- cts %>%
            filter(!!sym(cell_type_column) == cell_type) %>%
            filter(condition == condition_name) %>%
            filter(type == denoising) %>%
            filter(median_cts > expr_threshold) %>%
            pull(marker)
    } else {
        expr_mrk <- set_order
    }

    # Subset matrix
    m <- m[expr_mrk, expr_mrk]

    # Get the percentage of cells with proximity scores
    prox_pct <- prox_scores_summary %>%
        filter(!!sym(cell_type_column) == cell_type) %>%
        filter(condition == condition_name) %>%
        filter(type == denoising) %>%
        select(marker_1, marker_2, no_of_cells_with_prox_pct)
    prox_pct <- rbind(
        prox_pct,
        prox_pct %>% rename(marker_1 = marker_2, marker_2 = marker_1)
    )

    # Reorder markers
    if (is.null(set_order)) {
        hc <- hclust(dist(m), method = "ward.D2")
        marker_order <- hc$labels[hc$order]
        m <- m[marker_order, marker_order]
    } else {
        # Create a matrix with the size of set_order filled with zeros
        m0 <- matrix(0, nrow = length(set_order), ncol = length(set_order))
        rownames(m0) <- set_order
        colnames(m0) <- set_order

        # Find common markers between m and set_order
        common_markers <- intersect(rownames(m), set_order)

        # Fill m0 with values from m where markers exist in both matrices
        if (length(common_markers) > 0) {
            for (i in common_markers) {
                for (j in common_markers) {
                    m0[i, j] <- m[i, j]
                }
            }
        }

        # Use the ordered matrix for visualization
        m <- m0
    }

    # Limit values
    m[m > limit_values] <- limit_values
    m[m < -limit_values] <- -limit_values

    m_long <- m %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "marker_1") %>%
        pivot_longer(
            cols = -marker_1,
            names_to = "marker_2",
            values_to = "prox_score"
        )

    # Join with percentage of cells with proximity scores
    m_long <- m_long %>%
        left_join(prox_pct, by = c("marker_1", "marker_2"))

    # Order markers for plotting
    m_long$marker_1 <- factor(m_long$marker_1, levels = set_order)
    m_long$marker_2 <- factor(m_long$marker_2, levels = rev(set_order))

    heatmap_gradient <-
        colorRampPalette(
            c(
                "#1F395F",
                "#496389",
                "#728BB1",
                "#AABAD1",
                "#DFE5EE",
                "#FFFFFF",
                "#FFE0EA",
                "#E9AABF",
                "#CD6F8D",
                "#A23F5E",
                "#781534"
            )
        )

    p <- ggplot(
        m_long,
        aes(
            marker_1,
            marker_2,
            fill = prox_score,
            size = no_of_cells_with_prox_pct
        )
    ) +
        geom_point(shape = 21) +
        scale_size(range = c(0.1, 4), limits = c(0, 100)) +
        scale_x_discrete(position = "top") +
        # scale_fill_gradientn(
        #     colors = c(
        #         "#053061",
        #         "#2166AC",
        #         "#4393C3",
        #         "#92C5DE",
        #         "#D1E5F0",
        #         "#F7F7F7",
        #         "#FDDBC7",
        #         "#F4A582",
        #         "#D6604D",
        #         "#B2182B",
        #         "#67001F"
        #     ),
        #     limits = c(-limit_values, limit_values)
        # ) +
        scale_fill_gradientn(
            colors = heatmap_gradient(100),
            limits = c(-limit_values, limit_values)
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 0)) +
        coord_fixed() +
        labs(
            size = "Percentage of cells\n with proximity score",
            fill = paste0("Proximity Score\n", metric_name)
        )

    p <- p +
        labs(
            title = paste0(condition_name, ": ", cell_type),
            subtitle = paste0(
                "Only markers with median counts > ",
                expr_threshold,
                " are shown and heatmap is limited to ",
                limit_values,
                "\n",
                "denoising: ",
                denoising
            )
        )

    if (return_order) {
        return(marker_order)
    } else {
        return(p)
    }
}

create_metric_violins <- function(data, metrics) {
    # Calculate medians and max values for label positioning
    # Input: data frame and vector of metric column names
    # Returns: ggplot object with violin plots faceted by metric
    medians <- data[[]] %>%
        select(condition, all_of(metrics)) %>%
        pivot_longer(
            cols = all_of(metrics),
            names_to = "metric",
            values_to = "value"
        ) %>%
        group_by(metric, condition) %>%
        summarise(
            median_val = median(value),
            max_val = max(value),
            .groups = "drop"
        ) %>%
        group_by(metric) %>%
        mutate(label_pos = max(max_val) * 1.2) # Position labels 20% above max value

    # Create violin plots
    data[[]] %>%
        select(condition, all_of(metrics)) %>%
        pivot_longer(
            cols = all_of(metrics),
            names_to = "metric",
            values_to = "value"
        ) %>%
        # Create plot
        ggplot(aes(condition, value, color = condition, fill = condition)) +
        geom_violin(draw_quantiles = 0.5) +
        # Add median labels at consistent relative positions
        geom_text(
            data = medians,
            aes(y = label_pos, label = round(median_val, 1)),
            position = position_dodge(0.9),
            vjust = 0,
            color = "black"
        ) +
        facet_grid(metric ~ ., scales = "free") +
        scale_y_log10(expand = expansion(mult = c(0.05, 0.1))) + # Increased upper expansion to make room for labels
        theme_bw()
}


