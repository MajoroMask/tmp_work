
consite <- "/data1/suna/work/tmp_work/20230606_plot"
dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

library(rlang)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(scales)
library(fs)
library(ggstatsplot)

# func ----

my_breaks <- function(break_at) {
  break_fun <- function(limits) {
    inter <- limits[2] - limits[1]
    breaks <- c(limits[1] + inter * break_at)
    return(breaks)
  }
  return(break_fun)
}

# main ----

a <- new_environment()
a$pdi <- path("_i") %>% dir_create()
a$pdo <- path("_o") %>% dir_create()

## p1 ----

tb_p1_ori <-
  path(a$pdi, "20230601-CHO-IMBS-abundance normalized.csv") %>%
  vroom::vroom()
tb_p1 <-
  tb_p1_ori %>%
  transmute(
    MW = `MW [kDa]`,
    pI = `calc. pI`,
    sample_IMBS = `Abundance: F8: Sample, IMBS`,
    sample_STD = `Abundance: F10: Sample, STD`
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("sample_"),
    names_to = "group",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance)) %>%
  mutate(
    group = stringr::str_replace(group, "sample_", ""),
    log10_abundance = log10(abundance),
    log10_MW = log10(MW)
  )
p1 <-
  ggplot(tb_p1) +
  geom_point(
    aes(x = pI, y = log10(MW), size = abundance, fill = group, alpha = group),
    shape = 21, stroke = 0
  ) +
  scale_fill_manual(values = c("#f56025", "#3366cc")) +
  scale_size_continuous(
    range = c(2, 20),
    breaks = my_breaks(c(0.25, 0.5, 0.75)),
    labels = scales::label_number(scale_cut = cut_short_scale()),
    guide = guide_legend(override.aes = list(alpha = 1, stroke = 1))
  ) +
  scale_alpha_manual(values = c("IMBS" = 0.5, "STD" = 0.385)) +
  # facet_grid(rows = vars(group)) +
  theme_bw()
p1
ggsave(
  p1, filename = path(a$pdo, "20230606_bubble_cho_overlay.pdf"),
  width = 8, height = 4.5, scale = 1.5
)

p1_scatter <- ggstatsplot::ggscatterstats(
  tb_p1,
  x = pI,
  y = log10_MW,
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
  ylab = "Molecular weight (log10 scaled)"
)

p1_gs <- ggstatsplot::grouped_ggscatterstats(
  tb_p1,
  x = pI,
  y = log10_MW,
  grouping.var = group,
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
  ylab = "Molecular weight (log10 scaled)"
)
ggsave(
  p1_gs, file = path(a$pdo, "20230606_grouped_scatter_cho.pdf"),
  width = 16, height = 8
)

p1_ab <- ggstatsplot::ggbetweenstats(
  tb_p1,
  x = group,
  y = log10_abundance,
  bf.message = FALSE
)
p1_ab
ggsave(
  p1_ab, file = path(a$pdo, "20230606_box_cho.pdf"),
  height = 8, width = 6
)

## p1 new ----

tb_p1_new <-
  tb_p1_ori %>%
  transmute(
    MW = `MW [kDa]`,
    pI = `calc. pI`,
    abd_ratio = `Abundance Ratio: (IMBS) / (STD)`
  ) %>%
  filter(!is.na(abd_ratio)) %>%
  mutate(
    # group = stringr::str_replace(group, "sample_", ""),
    log2_abundance = log2(abd_ratio),
    log10_abundance = log10(abd_ratio),
    log10_MW = log10(MW)
  )
p1_new_n2pos <- ggstatsplot::gghistostats(
  # tb_p1_new %>% filter(log10_abundance != -2),
  tb_p1_new,
  x = log10_abundance,
  binwidth = 0.05
)
ggsave(
  p1_new_n2pos, filename = path(a$pdo, "20230606_ratio_histo.pdf"),
  width = 10, height = 4
)
p1_new_n2neg <- ggstatsplot::gghistostats(
  tb_p1_new %>% filter(log10_abundance != -2),
  x = log10_abundance,
  binwidth = 0.05
)
ggsave(
  p1_new_n2neg, filename = path(a$pdo, "20230606_ratio_histo_no0.01.pdf"),
  width = 10, height = 4
)

## p2 ----

tb_p2_ori <-
  path(a$pdi, "20230601-293T-IMBS-1ug-01.csv") %>%
  vroom::vroom()
tb_p2 <-
  tb_p2_ori %>%
  transmute(
    MW = `MW [kDa]`,
    pI = `calc. pI`,
    sample_IMBS = `Abundance: F5: Sample, IMBS`,
    sample_STD = `Abundance: F7: Sample, STD`,
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("sample_"),
    names_to = "group",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance)) %>%
  mutate(
    group = stringr::str_replace(group, "sample_", ""),
    log10_abundance = log10(abundance),
    log10_MW = log10(MW)
  )
p2 <-
  ggplot(tb_p2) +
  geom_point(
    aes(x = pI, y = log10(MW), size = abundance, fill = group, alpha = group),
    shape = 21, stroke = 0
  ) +
  scale_fill_manual(values = c("#f56025", "#3366cc")) +
  scale_size_continuous(
    range = c(2, 20),
    breaks = my_breaks(c(0.25, 0.5, 0.75)),
    labels = scales::label_number(scale_cut = cut_short_scale()),
    guide = guide_legend(override.aes = list(alpha = 1, stroke = 1))
  ) +
  scale_alpha_manual(values = c("IMBS" = 0.5, "STD" = 0.385)) +
  # facet_grid(rows = vars(group)) +
  theme_bw()
p2
ggsave(
  p2, filename = path(a$pdo, "20230606_bubble_293T_overlay.pdf"),
  width = 8, height = 4.5, scale = 1.5
)

p2_gs <- ggstatsplot::grouped_ggscatterstats(
  tb_p2,
  x = pI,
  y = log10_MW,
  grouping.var = group,
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
  ylab = "Molecular weight (log10 scaled)"
)
ggsave(
  p2_gs, file = path(a$pdo, "20230606_grouped_scatter_293T.pdf"),
  width = 16, height = 8
)

p2_ab <- ggstatsplot::ggbetweenstats(
  tb_p2,
  x = group,
  y = log10_abundance,
  bf.message = FALSE
)
p2_ab
ggsave(
  p2_ab, file = path(a$pdo, "20230606_box_293T.pdf"),
  height = 8, width = 6
)

## p3 ----

tb_p3_ori <-
  path(a$pdi, "20230330-ASK-ACvsDS-it100.csv") %>%
  vroom::vroom()
tb_p3 <-
  tb_p3_ori %>%
  transmute(
    MW = `MW [kDa]`,
    pI = `calc. pI`,
    sample_F151 = `Abundance: F151: Sample`,
    sample_F152 = `Abundance: F152: Sample`,
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("sample_"),
    names_to = "group",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance)) %>%
  mutate(
    group = stringr::str_replace(group, "sample_", ""),
    log10_abundance = log10(abundance),
    log10_MW = log10(MW)
  )
p3 <-
  ggplot(tb_p3) +
  geom_point(
    aes(x = pI, y = log10(MW), size = abundance, fill = group, alpha = group),
    shape = 21, stroke = 0
  ) +
  scale_fill_manual(values = c("#f56025", "#3366cc")) +
  scale_size_continuous(
    range = c(2, 20),
    breaks = my_breaks(c(0.25, 0.5, 0.75)),
    labels = scales::label_number(scale_cut = cut_short_scale()),
    guide = guide_legend(override.aes = list(alpha = 1, stroke = 1))
  ) +
  scale_alpha_manual(values = c("F151" = 0.75, "F152" = 0.585)) +
  # facet_grid(rows = vars(group)) +
  theme_bw()
p3
ggsave(
  p3, filename = path(a$pdo, "20230606_bubble_p3_overlay.pdf"),
  width = 8, height = 4.5, scale = 1.5
)

p3_gs <- ggstatsplot::grouped_ggscatterstats(
  tb_p3,
  x = pI,
  y = log10_MW,
  grouping.var = group,
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
  ylab = "Molecular weight (log10 scaled)"
)
ggsave(
  p3_gs, file = path(a$pdo, "20230606_grouped_scatter_p3.pdf"),
  width = 16, height = 8
)

p3_ab <- ggstatsplot::ggbetweenstats(
  tb_p3,
  x = group,
  y = log10_abundance,
  bf.message = FALSE
)
p3_ab
ggsave(
  p3_ab, file = path(a$pdo, "20230606_box_p3.pdf"),
  height = 8, width = 6
)
