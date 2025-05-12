heatmap_gradient <-
  colorRampPalette(
    c("#1F395F", "#496389", "#728BB1", "#AABAD1", "#DFE5EE", "#FFFFFF",
      "#FFE0EA", "#E9AABF", "#CD6F8D", "#A23F5E", "#781534"))

cherry_gradient <-
  colorRampPalette(
    c("gray90", "#FFE0EA", "#E9AABF", "#CD6F8D", "#A23F5E", "#781534"))

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
    "TCM CD4 T" = "#4A73C0",
    "TEM CD4 T" = "#224792",
    "Tregs" = "#1C3A76",
    "Naive CD8 T" = "#D0EDE6",
    "TCM CD8 T" = "#4AAF9D",
    "TEM CD8 T" = "#1B7E6F",
    "MAIT" = "#156559",
    "CD56dim NK" = "#A28EDB",
    "CD56bright NK" = "#866CCD",
    "Naive B" = "#F0DAC4",
    "Intermediate B" = "#C7A989",
    "Memory B" = "#917557",
    "mDC" = "#DA94C1",
    "pDC" = "#BB5391",
    "CD14 Mono" = "#E6BB43",
    "CD16 Mono" = "#AE8A1B",
    "Neutrophils" = "#CDCDCD",
    "Basophils" = "#797979",
    "Platelets" = "#7C2628",
    "gdT" = "#5A3E9E"
  )
