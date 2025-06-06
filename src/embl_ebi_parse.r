library(tidyverse)
library(here)

df <- read_csv(here("test.csv")) %>%
  rowwise() %>%
  mutate(type = case_when(
    feature_type == "SIGNAL" ~ "SP",
    grepl("NON CYTOPLASMIC.", description, fixed = TRUE) ~ "TM",
    TRUE ~ "OTHER"
  )) %>%
  filter(grepl("NON CYTOPLASMIC.", description, fixed = TRUE) |
    type == "SP")
