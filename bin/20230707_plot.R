
consite <- "/data1/suna/work/tmp_work/20230707_plot"
dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

library(rlang)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(scales)
library(fs)
library(ggstatsplot)

# utils ----

my_breaks <- function(break_at) {
  break_fun <- function(limits) {
    inter <- limits[2] - limits[1]
    breaks <- c(limits[1] + inter * break_at)
    return(breaks)
  }
  return(break_fun)
}

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
  path(a$pdi, "Ecoli-coverage.xlsx") %>%
  openxlsx::read.xlsx() %>%
  as_tibble()
tb_p <-
  tb_input %>%
  transmute(
    acc = Accession,
    mw = `MW.[kDa]`,
    log10_mw = log10(mw),
    pi = `calc..pI`,
    abundance = as.double(`Abundance.Ratio:.(F14)./.(F13)`),
    log2_abun = log2(abundance)
  )

## bubbles ----

p_bubble <-
  tb_p %>%
  filter(!is.na(log2_abun)) %>%
  arrange(abs(log2_abun)) %>%
  ggplot() +
  geom_point(
    aes(log10_mw, pi, color = log2_abun),
    size = 2.5
  ) +
  scale_color_gradient2() +
  scale_x_log10(n.breaks = 10) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
# p_bubble
ggsave(
  p_bubble, filename = path(a$pdo, "20230707_bubble.pdf"),
  width = 8, height = 6, scale = 1
)

cor(
  tb_p %>%
    filter(!is.na(log2_abun)) %>%
    select(pi, log2_abun) %>%
    as.matrix(),
  method = "pearson"
)

palette_hot <-
  c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000")
palette_cold <-
  c("#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d") %>% tail(-1)
p_bubble_pos <-
  tb_p %>%
  filter(!is.na(log2_abun)) %>%
  filter(log2_abun > 0) %>%
  arrange(abs(log2_abun)) %>%
  ggplot() +
  geom_point(
    aes(log10_mw, pi, color = log2_abun),
    size = 2.5
  ) +
  # scale_color_gradient() +
  # scale_color_viridis_c(option = "C") +
  scale_color_gradientn(colours = palette_hot) +
  scale_x_log10(n.breaks = 10) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
p_bubble_pos
p_bubble_neg <-
  tb_p %>%
  filter(!is.na(log2_abun)) %>%
  filter(log2_abun < 0) %>%
  arrange(abs(log2_abun)) %>%
  ggplot() +
  geom_point(
    aes(log10_mw, pi, color = log2_abun),
    size = 2.5
  ) +
  # scale_color_gradient() +
  # scale_color_viridis_c(option = "C") +
  scale_color_gradientn(colours = rev(palette_cold)) +
  scale_x_log10(n.breaks = 10) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
p_bubble_neg

ggsave(
  p_bubble_pos, filename = path(a$pdo, "20230707_bubble_pos.pdf"),
  width = 8, height = 6, scale = 1
)
ggsave(
  p_bubble_neg, filename = path(a$pdo, "20230707_bubble_neg.pdf"),
  width = 8, height = 6, scale = 1
)

## scatter ----

p_scatter <- ggstatsplot::ggscatterstats(
  tb_p,
  x = pi,
  y = log2_abun,
  # results.subtitle = FALSE,
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
  ylab = "Abundance Ratio F14/F13 (log2 scaled)"
)
ggsave(
  p_scatter,
  file = path(a$pdo, "20230707_scatter.pdf"),
  width = 10, height = 8
)

p_scatter_filtered <- ggstatsplot::ggscatterstats(
  tb_p %>% filter(abundance != 100, abundance != 0.01),
  x = pi,
  y = log2_abun,
  # results.subtitle = FALSE,
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
  ylab = "Abundance Ratio F14/F13 (log2 scaled)"
)
ggsave(
  p_scatter_filtered,
  file = path(a$pdo, "20230707_scatter_filtered.pdf"),
  width = 10, height = 8
)

p_scatter_filtered_pos <- ggstatsplot::ggscatterstats(
  data =
    tb_p %>%
    filter(abundance != 100, abundance != 0.01) %>%
    filter(abundance > 1),
  x = pi,
  y = log2_abun,
  # results.subtitle = FALSE,
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
  ylab = "Abundance Ratio F14/F13 (log2 scaled)"
)

p_scatter_filtered_neg <- ggstatsplot::ggscatterstats(
  data =
    tb_p %>%
    filter(abundance != 100, abundance != 0.01) %>%
    filter(abundance < 1),
  x = pi,
  y = log2_abun,
  # results.subtitle = FALSE,
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
  ylab = "Abundance Ratio F14/F13 (log2 scaled)"
)

## box, box ----

tb_pbox <-
  tb_input %>%
  transmute(
    acc = Accession,
    mw = `MW.[kDa]`,
    log10_mw = log10(mw),
    pi = `calc..pI`,
    abundance = as.double(`Abundance.Ratio:.(F14)./.(F13)`),
    log2_abun = log2(abundance),
    cat_abun = if_else(abundance > 1, "abundance > 1", "abundance < 1")
  )

p_box <- ggbetweenstats(
  data =
    tb_pbox %>%
    filter(!is.na(abundance)) %>%
    filter(abundance != 100, abundance != 0.01) %>%
    pipe_end(),
  x = cat_abun,
  y = pi,
  pairwise.comparisons = FALSE
)
ggsave(
  p_box,
  file = path(a$pdo, "20230707_box_filtered.pdf"),
  width = 10, height = 8
)
