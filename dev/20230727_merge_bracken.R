
consite <- "/data1/suna/work/tmp_work/20230727_merge_bracken"
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
a$pdi <- path("_i") %>% dir_create()
a$pdo <- path("_o") %>% dir_create()

# tb_bracken_1 <-
  path(a$pdi, "bracken_outputs_S_merge.txt") %>%
  vroom() %>%
  rowwise() %>%
  mutate(
    xit = sum(c_across(ends_with("_num")))
  ) %>%
  ungroup() %>%
  filter(xit != 0) %>%
  select(-xit) %>%
  pipe_end()
tb_bracken_2 <-
  path("/data1/suna/work/mag0724/output_mapping/Taxonomy/bracken/bracken_outputs_S.txt") %>%
  vroom()

sorted_colnames <-
  colnames(tb_bracken_1) %>%
  setdiff(c("name", "taxonomy_id", "taxonomy_lvl")) %>%
  sort() %>%
  {c(rbind(.[c(FALSE, TRUE)], .[c(TRUE, FALSE)]))} %>%
  c(
    c("name", "taxonomy_id", "taxonomy_lvl"),
    .,
    colnames(tb_bracken_2) %>% setdiff(c("name", "taxonomy_id", "taxonomy_lvl"))
  )

tb_merge <-
  full_join(
    tb_bracken_1, tb_bracken_2,
    by = c("name", "taxonomy_id", "taxonomy_lvl")
  ) %>%
  select(all_of(sorted_colnames)) %>%
  mutate(
    across(
      .cols = ends_with("_num"),
      .fns = ~ replace_na(.x, 0L)
    ),
    across(
      .cols = ends_with("_frac"),
      .fns = ~ replace_na(.x, 0)
    )
  )

# tb_test <-
#   tb_merge %>%
#   rowwise() %>%
#   mutate(
#     xit = sum(c_across(ends_with("_num")))
#   ) %>%
#   ungroup()
# tb_test %>% gsa(xit)

vroom::vroom_write(
  tb_merge,
  file = path(a$pdo, "20230727_bracken2_merged.txt"),
  delim = "\t"
)
