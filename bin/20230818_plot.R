
consite <- "/data1/suna/work/tmp_work/20230818_plot"
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
library(ggvenn)
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

tb_input <-
  path(a$pdi, "venn.xlsx") %>%
  openxlsx::read.xlsx() %>%
  as_tibble()
tb_p <-
  tb_input %>%
  rename_with(
    .cols = starts_with("Found.in.Sample:"),
    .fn = ~ str_replace(.x, ".*(\\[.*\\])\\.(.*)\\:.*", "\\1 \\2")
  ) %>%
  mutate(
    across(
      .cols = starts_with("["),
      .fns = ~ case_when(
        .x == "High" | .x == "Peak Found" ~ TRUE,
        .x == "Not Found" ~ FALSE,
        TRUE ~ NA
      )
    )
  )

p <-
  ggvenn::ggvenn(
    tb_p,
    columns = tb_p %>% select(starts_with("[")) %>% colnames(),
    digits = 2,
    fill_color = ggthemes::gdocs_pal()(4),
    stroke_color = "#333333"
  )
ggsave(
  p, filename = path(a$pdo, "venn.pdf"),
  width = 8, height = 4, scale = 2
)
