library(tidyverse)
library(here)

our_data_files <- list.files(
  "analysis/data/raw_data", pattern = ".csv$", full.names=TRUE)
#%>% str_subset(., "sun")

our_points <-
  map(our_data_files,
      read_csv)

names(our_points) <- our_data_files

our_points_tbl <-
  bind_rows(our_points,
            .id = "specimen")

our_points_tbl %>%
  filter(specimen %in% our_data_files[1:61]) %>%
  ggplot() +
  aes(V1, V2) +
  geom_point() +
  facet_wrap(~ specimen) +
  coord_fixed()
