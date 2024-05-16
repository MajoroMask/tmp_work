
consite <- "/data1/suna/work/tmp_work/20230706_plot"
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
  path(a$pdi, "MDCK-coverage-1ug.xlsx") %>%
  openxlsx::read.xlsx(sheet = "Sheet1") %>%
  as_tibble()
tb_p <-
  tb_input %>%
  transmute(
    acc = Accession,
    mw = `MW.[kDa]`,
    log10_mw = log10(mw),
    pi = `calc..pI`,
    F15 = `Found.in.Sample:.F15:.Sample` == "High",
    F16 = `Found.in.Sample:.F16:.Sample` == "High",
    F17 = `Found.in.Sample:.F17:.Sample` == "High",
    triple_threats = !(F15 & F16 & F17)
  ) %>%
  tidyr::pivot_longer(
    cols = F15:F17,
    names_to = "x",
    values_to = "inornot"
  ) %>%
  filter(inornot)

## original plots ----

p_mw <- ggstatsplot::ggbetweenstats(
  tb_p,
  x = x,
  y = log10_mw,
  xlab = "Peak found in",
  ylab = "Molecular weight (log10 scaled)",
  bf.message = FALSE,
  pairwise.comparisons = FALSE
)
p_pi <- ggstatsplot::ggbetweenstats(
  tb_p,
  x = x,
  y = pi,
  xlab = "Peak found in",
  ylab = "pI",
  bf.message = FALSE,
  pairwise.comparisons = FALSE
)

ggsave(p_mw, file = path(a$pdo, "p_mw.pdf"), height = 8, width = 8)
ggsave(p_pi, file = path(a$pdo, "p_pi.pdf"), height = 8, width = 8)

## triple threads annihilated ----

p_mw_tt <- ggstatsplot::ggbetweenstats(
  tb_p %>% filter(triple_threats),
  x = x,
  y = log10_mw,
  xlab = "Peak found in",
  ylab = "Molecular weight (log10 scaled)",
  bf.message = FALSE,
  pairwise.comparisons = FALSE
)
p_pi_tt <- ggstatsplot::ggbetweenstats(
  tb_p %>% filter(triple_threats),
  x = x,
  y = pi,
  xlab = "Peak found in",
  ylab = "pI",
  bf.message = FALSE,
  pairwise.comparisons = FALSE
)

ggsave(p_mw_tt, file = path(a$pdo, "p_mw_tt.pdf"), height = 8, width = 8)
ggsave(p_pi_tt, file = path(a$pdo, "p_pi_tt.pdf"), height = 8, width = 8)

## how about …… ----

# 完全不行，mw和abundance没有相关性

tb_p2 <-
  tb_input %>%
  transmute(
    acc = Accession,
    mw = `MW.[kDa]`,
    log10_mw = log10(mw),
    pi = `calc..pI`,
    ratio15 = as.double(`Abundance.Ratio:.(F15)./.(F17)`),
    ratio16 = as.double(`Abundance.Ratio:.(F16)./.(F17)`),
    log2_r15 = log2(ratio15),
    log2_r16 = log2(ratio16),
    F15 = `Found.in.Sample:.F15:.Sample` == "High",
    F16 = `Found.in.Sample:.F16:.Sample` == "High",
    F17 = `Found.in.Sample:.F17:.Sample` == "High",
    triple_threats = !(F15 & F16 & F17),
    enlist15 = (F15 & F17),
    enlist16 = (F16 & F17)
  )
cor(
  tb_p2 %>%
    filter(triple_threats) %>%
    select(log10_mw, log2_r15) %>%
    as.matrix(),
  use = "na.or.complete"
)
cor(
  tb_p2 %>%
    filter(triple_threats) %>%
    select(log10_mw, log2_r16) %>%
    as.matrix(),
  use = "na.or.complete"
)

tb_p2 %>%
  filter(enlist15) %>%
  arrange(abs(log2_r15)) %>%
  ggplot() +
  geom_point(
    # aes(log10_mw, pi, color = log2_r15),
    aes(mw, pi, color = log2_r15),
    size = 2.5
  ) +
  scale_color_gradient2() +
  scale_x_log10(n.breaks = 20) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
tb_p2 %>%
  filter(enlist16) %>%
  arrange(abs(log2_r16)) %>%
  ggplot() +
  geom_point(
    aes(log10_mw, pi, color = log2_r16),
    size = 2.5
  ) +
  scale_color_gradient2() +
  theme_bw()
