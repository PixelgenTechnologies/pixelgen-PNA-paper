
unique_annotations <-
  FetchData(pg_data_annotated, "cell_annotation") %>%
  pull(1) %>%
  unique() %>%
  sort()

annotated_cells %>% names %>% dput

cell_lineage <-
  list(Platelets = annotated_cells[["Platelets"]],
       Granulocytes = c(annotated_cells[["Basophils"]],
                        annotated_cells[["Neutrophils"]]),
       B = annotated_cells[["B"]],
       T = annotated_cells[["T"]],
       NK = annotated_cells[["NK"]],
       Mono = c(annotated_cells[["Classical Mono"]],
                annotated_cells[["Non-classical Mono"]],
                annotated_cells[["Intermediate Mono"]]),
       DC = c(annotated_cells[["pDC"]],
              annotated_cells[["cDC1"]],
              annotated_cells[["cDC2"]]))



sim_settings <-
  expand_grid(anno1 = names(cell_lineage),
              anno2 = names(cell_lineage)) %>%
  filter(anno1 < anno2) %>%
  # head() %>%
  unite(group_name, c(anno1, anno2), sep = "__", remove = FALSE) %>%
  group_by_all()

targeted_doublets <-
  sim_settings %>%
  group_split() %>%
  set_names(sim_settings$group_name) %>%
  pblapply(function(annos){
    annos <<-annos
    cells1 <- cell_lineage[[annos$anno1]]
    cells1 <- cells1[cells1 %in% colnames(pg_data_annotated)]
    # cell_annotation %>%
    # filter(cell_annotation == annos$anno1) %>%
    # pull(sample_component)

    cells2 <- cell_lineage[[annos$anno2]]
    cells2 <- cells2[cells2 %in% colnames(pg_data_annotated)]
    # cell_annotation %>%
    # filter(cell_annotation == annos$anno2) %>%
    # pull(sample_component)

    n_min <- pmin(length(cells1), length(cells2))
    sim_factor = 5

    sim_rate <- (n_min * sim_factor) / ncol(pg_data_annotated)


    PredictDoublets(pg_data_annotated,
                    ref_cells1 = cells1,
                    ref_cells2 = cells2,
                    npcs = 10,
                    simulation_rate = sim_rate,
                    n_neighbor = 50,
                    verbose = FALSE) %>%
      as_tibble(rownames = "sample_component")
  }) %>%
  bind_rows(.id = "group_name") %>%
  left_join(sim_settings)


targeted_doublets %>%
  group_by(anno1, anno2, doublet_prediction) %>%
  count() %>%
  group_by(anno1, anno2) %>%
  mutate(prc = n / sum(n) * 100) %>%
  filter(doublet_prediction == "doublet") %>%
  ggplot(aes(anno1, anno2, fill = prc)) +
  geom_tile() +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


plot_data <-
  pg_data_annotated %>%
  Embeddings("umap") %>%
  as_tibble(rownames = "sample_component") %>%
  left_join(FetchData(pg_data_annotated, c("seurat_clusters", "cell_annotation")) %>%
              as_tibble(rownames = "sample_component"),
            by = c("sample_component"))  %>%
  left_join(targeted_doublets) %>%
  arrange(desc(doublet_prediction))

plot_data %>%
  ggplot(aes(umap_1, umap_2, color = doublet_prediction)) +
  # geom_path() +
  geom_point(size = 2) +

  scale_color_manual(values = rev(yes_no_palette)) +
  coord_fixed() +
  facet_grid(anno1 ~ anno2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color = "Cell type",
       x = "UMAP 1",
       y = "UMAP 2")



plot_data %>%
  ggplot(aes(umap_1, umap_2, color = -log10(doublet_p_adj) )) +
  # geom_path() +
  geom_point(size = 2) +
  scale_color_gradientn(colors = cherry_gradient(100)) +
  coord_fixed() +
  facet_grid(anno1 ~ anno2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color = "Cell type",
       x = "UMAP 1",
       y = "UMAP 2")
