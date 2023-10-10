
consite <- "/data1/suna/work/tmp_work/20231010_plot"
dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

library(rlang)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(scales)
library(fs)
library(ggstatsplot)
library(ComplexUpset)

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

# main ----

a <- new_environment()
a$pdi <- path("_i") %>% dir_create()
a$pdo <- path("_o") %>% dir_create()

tb_input <-
  path(a$pdi, "mw_pi_formatted.xlsx") %>%
  openxlsx::read.xlsx() %>%
  as_tibble()
tb_p <-
  tb_input %>%
  transmute(
    acc = Accession,
    sample = sample,
    mw = `MW.[kDa]`,
    log2_mw = log2(mw),
    log10_mw = log10(mw),
    pi = `calc..pI`
  )

## upset ----

tb_grouping <-
  tb_p %>%
  mutate(dummy_fill = TRUE) %>%
  tidyr::pivot_wider(
    id_cols = acc,
    names_from = sample,
    values_from = dummy_fill,
    values_fill = FALSE
  )
samples <- colnames(tb_grouping) %>% tail(-1)

p_upset <- ComplexUpset::upset(
  tb_grouping,
  intersect = samples,
  min_size = 4
)

ggsave(
  p_upset, filename = path(a$pdo, "20231010_upset.pdf"),
  width = 12, height = 8
)

## scatter ----

l_scatters <-
  tb_p %>%
  group_by(sample) %>%
  group_map(
    .f = ~
      ggstatsplot::ggscatterstats(
        .x,
        x = pi,
        y = log2_mw,
        results.subtitle = FALSE,
        xsidehistogram = list(
          fill = "#F28E2B",
          color = "black", na.rm = TRUE, bins = 50
        ),
        ysidehistogram = list(
          fill = "#4E79A7",
          color = "black", na.rm = TRUE, bins = 50
        ),
        point.args = list(size = 3, alpha = 0.4, stroke = 0, color = "#E15759"),
        smooth.line.args = list(se = FALSE),
        xlab = "pI",
        ylab = "molecular weight (log2)"
      )
  )
purrr::walk2(
  .x = l_scatters,
  .y = tb_p %>% group_by(sample) %>% group_keys() %>% pull(1),
  .f = ~ ggsave(
    .x,
    file = path(a$pdo, glue::glue("20231010_scatter_{.y}.pdf")),
    width = 10, height = 8
  )
)

tb_scatter_union <-
  tb_p %>%
  group_by(acc) %>%
  slice_head(n = 1L) %>%
  ungroup()

p_scatter_union <- ggstatsplot::ggscatterstats(
  tb_scatter_union,
  x = pi,
  y = log2_mw,
  results.subtitle = FALSE,
  xsidehistogram = list(
    fill = "#F28E2B",
    color = "black", na.rm = TRUE, bins = 50
  ),
  ysidehistogram = list(
    fill = "#4E79A7",
    color = "black", na.rm = TRUE, bins = 50
  ),
  point.args = list(size = 3, alpha = 0.4, stroke = 0, color = "#E15759"),
  smooth.line.args = list(se = FALSE),
  xlab = "pI",
  ylab = "molecular weight (log2)"
)
ggsave(
  p_scatter_union,
  file = path(a$pdo, "20231010_scatter_union.pdf"),
  width = 10, height = 8
)
