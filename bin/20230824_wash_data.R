
consite <- "/data1/suna/work/mag0818"
R_proj <- "~/proj/tmp_work"

dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

renv::activate(R_proj)

library(vroom)
library(rlang)
library(purrr)
library(stringr)
library(tidyr)
library(forcats)
library(ggplot2)
library(ggstatsplot)
library(ggpubr)
library(ggthemes)
library(dplyr)
library(fs)

# utils ----

gsa <- function(.data, ..., .sort = TRUE) {
  requireNamespace("dplyr", quietly = TRUE)
  .data %>%
    group_by(...) %>%
    tally(sort = .sort) %>%
    ungroup() %>%
    # mutate(pct = round(n / sum(n), digits = 3))
    mutate(pct = scales::percent(round(n / sum(n), digits = 3)))
}

pipe_end <- function(x) x

# IO ----

a <- new_environment()
a$pdi <- path("output_dump_hash/Taxonomy/virus_mapping") %>% path_real()
a$pdo <- a$pdi

tb_xits_from_wyl <-
  path(consite, "Cluster_rep_taxid_taxonomy_check_host_Gtype.xls") %>%
  vroom::vroom()

# main ----

tb_covstat_input <-
  path(a$pdi, "bbmap_covstats.csv") %>%
  vroom::vroom(show_col_types = FALSE)
