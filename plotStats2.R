library(tidyverse)
args <- commandArgs(TRUE)
readFile <- function(fileName) {
    fileName |>
    readr::read_tsv(col_names=c("ctg", "start", "end", "mol_cov", "middle_mol_cov", "mol_len", "read_dens", "# start", "# end")) |>
    dplyr::select(-c("ctg", "start", "end")) |>
    tidyr::pivot_longer(cols = everything(), names_to = "data", values_to = "values") |>
    dplyr::mutate(data = factor(data))
}
dplyr::bind_cols(readFile(args[1]),
                 readFile(args[2]) |> dplyr::rename(dataScore = data, valuesScore = values)) |>
   ggplot(aes(x = values, y = valuesScore)) +
       #geom_bin2d() +
       geom_point() +
       facet_wrap(vars(data), scales = "free") + xlim(0, 1000) + ylim(0, 0.1)
