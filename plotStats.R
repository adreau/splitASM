library(tidyverse)
readr::read_tsv("tmp.out", col_names=c("ctg", "start", "end", "mol_cov", "middle_mol_cov", "mol_len", "read_dens", "# start", "# end")) |>
  dplyr::select(-c("ctg", "start", "end")) |>
  tidyr::pivot_longer(cols = everything(), names_to = "data", values_to = "values") |>
  dplyr::mutate(data = factor(data)) |>
  ggplot(aes(x = values)) +
    geom_histogram() +
    facet_wrap(vars(data), scales = "free") + xlim(0, 0.1) #  + xlim(0, 250)

