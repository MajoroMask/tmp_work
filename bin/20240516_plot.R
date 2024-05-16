
# 画散点图 ----

## init ----

a <- new.env(parent = emptyenv())
a$path_project <- "~/proj/tmp_work/"
a$pwd <- "/data1/suna/work/tmp_work/20240516_plot"
suppressWarnings(dir.create(a$pwd))
setwd(a$pwd)

library(dplyr)
library(purrr)
library(stringr)
library(fs)
library(ggstatsplot)

a$pdi <- path(a$pwd, "_i") |> dir_create()
a$pdo <- path(a$pwd, "_o") |> dir_create()

# func ----

gsa <- function(.data, ..., .sort = TRUE) {
  requireNamespace("dplyr", quietly = TRUE)
  .data |>
    group_by(...) |>
    tally(sort = .sort) |>
    ungroup() |>
    # mutate(pct = round(n / sum(n), digits = 3))
    mutate(pct = scales::percent(round(n / sum(n), digits = 3)))
}

## IO ----

tb_input_all <-
  map2_dfr(
    .x = fs::dir_ls(a$pdi, regexp = "xlsx$"),
    .y = fs::dir_ls(a$pdi, regexp = "xlsx$") |>
      path_file() |>
      str_replace(".*-(Sample_\\d+)-.*", "\\1") |>
      str_to_lower(),
    .f = ~
      openxlsx::read.xlsx(.x, sep.names = " ") |>
      as_tibble() |>
      mutate(sample = .y, .before = 1)
  )

tb_gsa <-
  tb_input_all |>
  filter(`Marked as` == "CHO") |>
  gsa(Accession) |>
  mutate() |>
  vroom::vroom_write(
    file = path(a$pdo, "protein_frequency.csv"),
    delim = ","
  )

## plot ----

tb_p <-
  tb_input_all |>
  filter(`Marked as` == "CHO") |>
  group_by(Accession) |>
  summarise(
    n = n(),
    pct = scales::percent(round(n / 26L, digits = 3)),
    pi = `calc. pI`[1],
    mw = log10(`MW [kDa]`[1]),
    .groups = "drop"
  ) |>
  ungroup() |>
  arrange(n) |>
  mutate(acc = fct_inorder(Accession)) |>
  mutate(
    color =
      RColorBrewer::brewer.pal(n = 11, "Spectral") |>
      _[c(10, 6, 2)] |>
      colorRampPalette() %>%
      {.(26)} %>%
      `[`(n)
  )
vroom::vroom_write(
  tb_p |> select(-acc, -color),
  file = path(a$pdo, "protein_frequency.csv"),
  delim = ","
)

p_bw <-
  ggscatterstats(
    data = tb_p,
    x = pi,
    y = mw,
    xlab = "pI",
    ylab = "MW (log10)",
    type = "np",
    pairwise.comparisons = FALSE,
    results.subtitle = FALSE,
    package = "ggthemes",
    palette = "gdoc",
    point.args = list(color = "black", size = 3, alpha = 0.4, stroke = 0),
    xsidehistogram.args = list(
      fill = "#009E73", color = "black",
      na.rm = TRUE, bins = 50
    ),
    ysidehistogram.args = list(
      fill = "#D55E00", color = "black",
      na.rm = TRUE, bins = 50
    ),
    smooth.line.args = list(
      linewidth = 1.5,
      color = "blue",
      method = "glm",
      formula = y ~ x
    )
  )
p_colored <-
  ggscatterstats(
    data = tb_p,
    x = pi,
    y = mw,
    xlab = "pI",
    ylab = "MW (log10)",
    type = "np",
    pairwise.comparisons = FALSE,
    results.subtitle = FALSE,
    package = "ggthemes",
    palette = "gdoc",
    point.args = list(color = tb_p$color, size = 3, alpha = 1, stroke = 0),
    xsidehistogram.args = list(
      fill = "#009E73", color = "black",
      na.rm = TRUE, bins = 50
    ),
    ysidehistogram.args = list(
      fill = "#D55E00", color = "black",
      na.rm = TRUE, bins = 50
    ),
    smooth.line.args = list(
      linewidth = 1.5,
      color = "blue",
      method = "glm",
      formula = y ~ x
    )
  )

walk2(
  .x = list(p_bw, p_colored),
  .y = c("p_bw", "p_colored"),
  .f = ~
    ggsave(
      .x, filename = path(a$pdo, glue("{.y}.pdf")),
      width = 8, height = 8, dpi = 300
    )
)
