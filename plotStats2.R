library(tidyverse)
args <- commandArgs(TRUE)
readFile1 <- function(fileName) {
    fileName |>
    readr::read_tsv(col_names=c("ctg", "start", "end", "mol_cov", "mol_len", "read_dens", "# start", "# end")) |>
    dplyr::select(-c("ctg", "start", "end")) |>
    tidyr::pivot_longer(cols = everything(), names_to = "data", values_to = "values") |>
    dplyr::mutate(data = factor(data))
}
readFile2 <- function(fileName) {
    fileName |>
    readr::read_tsv(col_names=c("ctg", "start", "end", "score")) |>
    dplyr::pull(score) |>
    rep(each = 5)
}

args[1] |>
  readFile1() |>
  dplyr::mutate(score = readFile2(args[2])) |>
   ggplot(aes(x = values, y = score)) +
       #geom_bin2d() +
       geom_point() +
       facet_wrap(vars(data), scales = "free")
