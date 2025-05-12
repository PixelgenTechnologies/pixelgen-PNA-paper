
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
