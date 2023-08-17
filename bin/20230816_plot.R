
consite <- "/data1/suna/work/tmp_work/20230816_plot"
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

# IO ----

a <- new_environment()
a$pdi <- path("_i") %>% dir_create()
a$pdo <- path("_o") %>% dir_create()

tb_input <-
  path(a$pdi, "附表1.HCCF蛋白鉴定及定量列表(1).xlsx") %>%
  openxlsx::read.xlsx() %>%
  as_tibble()

# main ----

tb_p <-
  tb_input %>%
  select(
    acc = Accession,
    abun1 = `DS2-HCCF.abuandance`,
    abun2 = `DS3-HCCF.abundance`,
    ratio = `Ratio(DS2/DS3)`
  ) %>%
  filter(
    !is.na(abun1),
    !is.na(abun2)
  ) %>%
  mutate(
    color =
      case_when(
        log2(ratio) > 2 ~ "up",
        log2(ratio) < -2 ~ "down",
        TRUE ~ "not sign"
      ) %>%
      factor(levels = c("not sign", "down", "up"))
  ) %>%
  arrange(color) %>%
  mutate(color = factor(color, levels = c("up", "not sign", "down")))

major_breaks_at <- c(0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000)
minor_breaks_at <- c(25, 75, 250, 750, 2500, 7500, 25000, 75000)
gather_breaks <- c(major_breaks_at, minor_breaks_at) %>% sort()

p <-
  ggplot(tb_p) +
  geom_abline(slope = 1, intercept = 0, color = alpha("#aaaaaa", 0.3)) +
  geom_abline(slope = 1, intercept = 2, color = alpha("#d7191c", 0.3)) +
  geom_abline(slope = 1, intercept = -2, color = alpha("#2c7bb6", 0.3)) +
  geom_point(
    aes(abun2, abun1, color = color, size = color),
    shape = 19, stroke = 0
  ) +
  # scale_x_continuous(
  #   trans = "log2",
  #   breaks = major_breaks_at,
  #   minor_breaks = minor_breaks_at
  # ) +
  # scale_y_continuous(
  #   trans = "log2",
  #   breaks = major_breaks_at,
  #   minor_breaks = minor_breaks_at
  # ) +
  scale_x_continuous(
    trans = "log2",
    breaks = gather_breaks
  ) +
  scale_y_continuous(
    trans = "log2",
    breaks = gather_breaks
  ) +
  scale_color_manual(
    values = c(
      "up" = "#d7191c",
      "not sign" = alpha("#aaaaaa", 1),
      "down" = "#2c7bb6"
    ),
    guide = guide_legend(
      title = "Ratio(DS2/DS3)",
      override.aes = list(size = 2)
    )
  ) +
  scale_size_manual(
    values = c("up" = 2, "not sign" = 1.5, "down" = 2),
    guide = "none"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(hjust = 1, angle = 45)
  ) +
  coord_fixed() +
  labs(
    x = "DS3-HCCF abuandance (ng, log2 scaled)",
    y = "DS2-HCCF abuandance (ng, log2 scaled)"
  )
ggsave(
  p, filename = path(a$pdo, "abun_and_ratio.pdf"),
  height = 6, width = 8
)
