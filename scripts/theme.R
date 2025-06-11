heatmap_gradient <-
  colorRampPalette(
    c("#1F395F", "#496389", "#728BB1", "#AABAD1", "#DFE5EE", "#FFFFFF",
      "#FFE0EA", "#E9AABF", "#CD6F8D", "#A23F5E", "#781534"))

cherry_gradient <-
  colorRampPalette(
    c("gray90", "#FFE0EA", "#E9AABF", "#CD6F8D", "#A23F5E", "#781534"))


cherry_light_gradient <-
  colorRampPalette(
    c("white", "#FFE0EA", "#E9AABF", "#CD6F8D", "#A23F5E", "#781534"))

yes_no_palette <-
  c("#DAD6D7", "#C86584")


theme_default <-
  function(base_size = 10) {
    theme_bw(base_size = base_size) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
  }



pbmc_cell_type_palette <-
  c(
    "Naive CD4 T" = "#B9CDED",
    "CD4 TSCM" = "#92B0E0",
    "CD4 TCM" = "#6D92D1",
    "CD4 TEM" = "#4A73C0",
    "CD4 TEFF" = "#224792",
    "CD4 TEMRA" = "#1C3A76",
    "Treg" = "#2955AE",
    "Naive CD8 T" = "#D0EDE6",
    "CD8 TSCM" = "#A2DACE",
    "CD8 TCM" = "#75C5B5",
    "CD8 TEM" = "#4AAF9D",
    "CD8 TEFF" = "#209785",
    "CD8 TEMRA" = "#1B7E6F",
    "MAIT" = "#156559",
    "DPT" = "#1B7E9F",
    "DNT" = "#7B787F",
    "CD56dim NK" = "#A28EDB",
    "CD56bright NK" = "#866CCD",
    "NKT" = "#6C4ABD",
    "Naive B" = "#F0DAC4",
    "Memory B" = "#917557",
    "Plasma cells" = "#DE9982",
    "cDC1" = "#DA94C1",
    "cDC2" = "#BB5391",
    "pDC" = "#9C4579",
    "Classical Mono" = "#F0C966",
    "Intermediate Mono" = "#DAAC22",
    "Non-classical Mono" = "#836714",
    "Neutrophils" = "#CDCDCD",
    "Basophils" = "#797979",
    "Platelets" = "#7C2628",
    "gdT" = "#5A3E9E"
  )

