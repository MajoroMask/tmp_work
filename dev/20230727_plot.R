
consite <- "/data1/suna/work/tmp_work/20230727_protein_cov"
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

# function ----

add_pep_pos_info <- function(tb) {
  tb_out <-
    tb %>%
    mutate(
      pep_pos =
        str_replace(Identification, ".*([A-Z]\\d+-.*[A-Z]\\d+) = .*", "\\1"),
      pep_start =
        str_replace(pep_pos, "[A-Z](\\d+)-[A-Z](\\d+)", "\\1") %>%
        as.integer(),
      pep_end =
        str_replace(pep_pos, "[A-Z](\\d+)-[A-Z](\\d+)", "\\2") %>%
        as.integer()
    )
  return(tb_out)
}

calculate_covstat <- function(tb) {
  tb_out <-
    tb %>%
    group_by(name, .drop = FALSE) %>%
    summarise(
      length = length(cov),
      covered_aa = sum(cov != 0),
      covered_pct = covered_aa / length,
      mean_fold = mean(cov),
      median_fold = median(cov),
      sd_fold = sd(cov),
      covered_mean_fold = cov[cov != 0L] %>% mean(),
      covered_median_fold = cov[cov != 0L] %>% median(),
      covered_sd_fold = cov[cov != 0L] %>% sd()
    ) %>%
    ungroup()
  return(tb_out)
}

# IO ----

a <- new_environment()
a$pdi <- path("_i") %>% dir_create()
a$pdo <- path("_o") %>% dir_create()

ori_fa <-
  path(a$pdi, "GOI4.fasta") %>%
  readr::read_lines()
tb_a <-
  path(a$pdi, "YMJN_Chy.csv") %>%
  vroom() %>%
  add_pep_pos_info()
tb_b <-
  path(a$pdi, "YMJN_Typsin_filter.csv") %>%
  vroom() %>%
  add_pep_pos_info()

# main ----

## data washing ----

tb_fa <-
  map2_dfr(
    .x =
      ori_fa[c(TRUE, FALSE)] %>%
      str_replace("^>", "") %>%
      str_replace("VSV-G", "VSVG"),
    .y = ori_fa[c(FALSE, TRUE)],
    .f = function(name, seq) {
      tb_out <-
        tibble(
          name = name,
          aa = str_split(seq, "") %>% unlist()
        ) %>%
        mutate(pos = row_number())
      return(tb_out)
    }
  ) %>%
  mutate(
    cov_Chy = 0,
    cov_Trypsin = 0
  )
tb_cov_flat <-
  tb_fa %>%
  reduce2(
    .init = .,
    .x = tb_a$pep_pos,
    .y = tb_a$Protein,
    .f = function(tb, pep_pos, protein) {
      start <-
        str_replace(pep_pos, "[A-Z](\\d+)-[A-Z](\\d+)", "\\1") %>%
        as.integer()
      end <-
        str_replace(pep_pos, "[A-Z](\\d+)-[A-Z](\\d+)", "\\2") %>%
        as.integer()
      tb_out <-
        tb %>%
        mutate(
          cov_Chy = if_else(
            name == protein & pos >= start & pos <= end,
            cov_Chy + 1L,
            cov_Chy
          )
        )
      return(tb_out)
    }
  ) %>%
  reduce2(
    .init = .,
    .x = tb_b$pep_pos,
    .y = tb_b$Protein,
    .f = function(tb, pep_pos, protein) {
      start <-
        str_replace(pep_pos, "[A-Z](\\d+)-[A-Z](\\d+)", "\\1") %>%
        as.integer()
      end <-
        str_replace(pep_pos, "[A-Z](\\d+)-[A-Z](\\d+)", "\\2") %>%
        as.integer()
      tb_out <-
        tb %>%
        mutate(
          cov_Trypsin = if_else(
            name == protein & pos >= start & pos <= end,
            cov_Trypsin + 1L,
            cov_Trypsin
          )
        )
      return(tb_out)
    }
  )

tb_stat_chy <-
  tb_cov_flat %>%
  mutate(
    name = factor(name) %>% fct_inorder(),
    cov = cov_Chy
  ) %>%
  calculate_covstat()
tb_stat_Trypsin <-
  tb_cov_flat %>%
  mutate(
    name = factor(name) %>% fct_inorder(),
    cov = cov_Trypsin
  ) %>%
  calculate_covstat()
tb_stat_combined <-
  tb_cov_flat %>%
  mutate(
    name = factor(name) %>% fct_inorder(),
    cov = cov_Trypsin + cov_Chy
  ) %>%
  calculate_covstat()

## output ----

wb <- openxlsx::createWorkbook(creator = "shenkebio")
openxlsx::addWorksheet(wb, sheetName = "Chy")
openxlsx::writeData(wb, sheet = "Chy", x = tb_stat_chy)
openxlsx::addWorksheet(wb, sheetName = "Trypsin")
openxlsx::writeData(wb, sheet = "Trypsin", x = tb_stat_Trypsin)
openxlsx::addWorksheet(wb, sheetName = "combined")
openxlsx::writeData(wb, sheet = "combined", x = tb_stat_combined)
openxlsx::saveWorkbook(wb, file = path(a$pdo, "covstat.xlsx"), overwrite = TRUE)

## cov plot ----

tb_p <-
  tb_cov_flat %>%
  mutate(
    cov_combined = cov_Chy + cov_Trypsin
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("cov_"),
    names_to = "group",
    values_to = "cov",

  ) %>%
  mutate(
    name = factor(name) %>% fct_inorder(),
    group =
      str_c(str_replace(group, "^cov_", ""), " - ", name, sep = "") %>%
      factor() %>%
      fct_inorder(),
    size = if_else(cov == 0, "zero", "non-zero")
  )
p <-
  ggplot(tb_p, aes(pos, cov, group = name, color = name)) +
  geom_line() +
  geom_point(aes(size = size)) +
  scale_color_tableau(guide = "none") +
  scale_size_manual(
    values = c("zero" = 0.25, "non-zero" = 1.5),
    guide = "none"
  ) +
  facet_wrap(
    facets = vars(group),
    ncol = 3,
    scales = "free"
  ) +
  theme_bw() +
  labs(
    x = "AA position",
    y = "Depth"
  )
ggsave(
  p, filename = path(a$pdo, "depth.pdf"),
  height = 8, width = 18, scale = 1.25
)
