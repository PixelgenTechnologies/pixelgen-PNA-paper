
library(jsonlite)
library(tidyverse)

files <- list.files("data", pattern = "json")

# read a json file
json_data <-
  files %>%
  set_names() %>%
  lapply(function(x) jsonlite::fromJSON(paste0("data/", x)))

plot_data <-
  json_data %>%
  map(. %>%
        {.$pre_filtering_component_sizes} %>%
        unlist() %>%
        enframe("n_umi", "number_of_cells")) %>%
  bind_rows(.id = "file") %>%
  mutate(n_umi = as.numeric(n_umi))


plot_data %>%
  arrange(-n_umi) %>%
  group_by(file) %>%
  mutate(cumsum = cumsum(number_of_cells)) %>%
  ggplot(aes(x = cumsum, y = n_umi, color = file)) +
  geom_step(orientation = "y") +
  geom_hline(yintercept = 8000, linetype = "dashed", color = "red") +
  # geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Cumulative Number of Cells", y = "Number of UMIs", color = "Sample")

ggsave("results/pbmc_results/umi_distribution_before_filtering.pdf", width = 6, height = 3.5)
