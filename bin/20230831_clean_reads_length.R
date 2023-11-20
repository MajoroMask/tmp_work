
consite <- "/data1/suna/work/tmp_work/20230831_clean_reads_length"
R_proj <- "~/proj/tmp_work"

dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

renv::activate(R_proj)

library(pheatmap)
library(vroom)
library(rlang)
library(purrr)
library(glue)
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

# func ----

read_overall_stats <- function(smp, path_input) {
  tb <-
    path_input %>%
    vroom(
      delim = "\t", n_max = 8, col_names = c("x", "y"),
      show_col_types = FALSE
    ) %>%
    mutate(
      x = str_replace(x, "^\\#(.*)\\:$", "\\1") %>% str_to_lower()
    ) %>%
    pivot_wider(names_from = x, values_from = y) %>%
    mutate(sample = smp, .before = 1)
  return(tb)
}

read_detail_stats <- function(smp, path_input) {
  tb <-
    path_input %>%
    vroom(
      delim = "\t", skip = 9,
      show_col_types = FALSE
    ) %>%
    rename(
      read_length = `#Length`,
      n_reads = reads,
      n_bases = bases,
    ) %>%
    mutate(
      across(
        .cols = matches("pct_"),
        .fns = ~
          .x %>%
          str_replace("\\%$", "") %>%
          as.double() %>%
          `/`(100)
      )
    ) %>%
    mutate(sample = smp, .before = 1)
  return(tb)
}

# IO ----

a <- new_environment()
a$pdo <- path("_o") %>% dir_create()

tb_samples_ori <-
  path("samplesheet.csv") %>%
  vroom(show_col_types = FALSE)
tb_samples <-
  tb_samples_ori %>%
  mutate(
    type = case_when(
      str_detect(sample, "NTC_DNA") ~ "NTC_DNA",
      str_detect(sample, "NTC") ~ "NTC_RNA",
      str_detect(sample, "293T_RNA_BVDV") ~ "293T_RNA_BVDV",
      str_detect(sample, "293T_DNA") ~ "293T_DNA",
      str_detect(sample, "293T_RNA") ~ "293T_RNA",
      # str_detect(sample, "MDBK_DNA") ~ "MDBK_DNA",
      # str_detect(sample, "MDBK_RNA") ~ "MDBK_RNA",
      str_detect(sample, "MDBK") ~ "MDBK",
      str_detect(sample, "BVDV") ~ "BVDV",
      TRUE ~ NA_character_
    ),
    is_ntc = str_detect(sample, "NTC"),
    path_rlstat = path("rl_output", glue("{sample}.rl.txt")) %>% path_real()
  ) %>%
  arrange(type, is_ntc, sample) %>%
  mutate(sample = factor(sample) %>% fct_inorder())

tb_rlstat <-
  purrr::map2_dfr(
    .x = tb_samples$sample,
    .y = tb_samples$path_rlstat,
    .f = read_overall_stats
  )
tb_rl_detail <-
  purrr::map2_dfr(
    .x = tb_samples$sample,
    .y = tb_samples$path_rlstat,
    .f = read_detail_stats
  )

tb_samples_extended <-
  tb_samples %>%
  left_join(tb_rlstat, by = "sample") %>%
  mutate(sample = factor(sample) %>% fct_inorder())

## overall distribution ----

tb_rl_boxstat <-
  tb_rl_detail %>%
  group_by(sample) %>%
  group_map(
    .f = function(tb, smp) {
      x <- rep(tb$read_length, times = tb$n_reads)
      uni_x <- tb$read_length
      qtl <- quantile(x, probs = seq(0, 1, 0.1))
      quad_qtl <- quantile(x, c(0.25, 0.75))
      iqr <- quad_qtl[2] - quad_qtl[1]
      lower_whisker <- min(uni_x[uni_x >= quad_qtl[1] - 1.5 * iqr])
      upper_whisker <- max(uni_x[uni_x <= quad_qtl[2] + 1.5 * iqr])
      tb_out <-
        tibble(
          name = sprintf("Q%02d", seq(0, 100, 10)),
          qtl = qtl
        ) %>%
        pivot_wider(names_from = name, values_from = qtl) %>%
        mutate(
          sample = tb$sample[1],
          q1 = quad_qtl[1],
          q3 = quad_qtl[2],
          lower = lower_whisker,
          upper = upper_whisker
        )
      return(tb_out)
    },
    .keep = TRUE
  ) %>%
  bind_rows()

tb_p_overall_dist <-
  tb_samples_extended %>%
  select(sample, type, is_ntc, median, mean = avg, min, max) %>%
  left_join(tb_rl_boxstat, by = "sample") %>%
  mutate(sample = factor(sample) %>% fct_inorder())

p_overall_dist <-
  ggplot(tb_p_overall_dist) +
  geom_boxplot(
    aes(
      x = sample,
      lower = q1,
      upper = q3,
      middle = median,
      ymin = lower,
      ymax = upper,
      color = is_ntc,
      fill = type
    ),
    stat = "identity",
    alpha = I(0.5)
  ) +
  scale_y_continuous(n.breaks = 10L) +
  scale_color_manual(values = c("#d53e4f", "#3288bd")) +
  scale_fill_gdocs() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 20, l = 100)
  ) +
  labs(y = "read length")
ggsave(
  p_overall_dist, filename = path(a$pdo, "overall_dist.svg"),
  width = 12, height = 8
)

## accumulated percentage ----

tb_p_cum_pct <-
  tb_rl_detail %>%
  left_join(
    tb_samples_extended %>% select(sample, type, is_ntc),
    by = "sample"
  ) %>%
  mutate(sample = factor(sample, levels = levels(tb_samples_extended$sample)))

l_p <-
  tb_p_cum_pct %>%
  group_by(type) %>%
  group_map(
    .f = ~
      ggplot(.x) +
      # geom_point(
      #   aes(read_length, y = cum_pct_reads, color = sample),
      #   size = 1
      # ) +
      geom_line(
        aes(read_length, y = cum_pct_reads, color = sample),
        alpha = I(0.75), linewidth = 1
      ) +
      scale_x_reverse(breaks = seq(0, 150, 10)) +
      scale_y_continuous(
        limits = c(0, 1),
        n.breaks = 5L,
        labels = scales::percent
      ) +
      scale_color_tableau(palette = "Tableau 10") +
      facet_wrap(facets = vars(type), nrow = 2L) +
      theme_bw() +
      theme(
        # legend.position = "bottom",
        strip.background = element_rect(fill = alpha("black", 0))
      ) +
      labs(y = "Accumulated read percentage"),
    .keep = TRUE
  )
purrr::walk2(
  .x = l_p,
  .y =
    tb_p_cum_pct %>%
    group_by(type) %>%
    dplyr::group_keys() %>%
    pull(type),
  .f = ~ ggsave(
    .x, filename = path(a$pdo, glue("accu_pct_read_{.y}.svg")),
    width = 10, height = 6
  )
)

## heatmap ----

m_heatmap <-
  tb_p_overall_dist %>%
  select(sample, starts_with("Q", ignore.case = FALSE)) %>%
  tibble::column_to_rownames(var = "sample") %>%
  as.matrix()
svg(file = path(a$pdo, "heatmap.svg"), width = 16, height = 8)
pheatmap::pheatmap(
  m_heatmap,
  # border_color = NA,
  cluster_cols = FALSE,
  annotation_row =
    tb_p_overall_dist %>%
    select(sample, type) %>%
    tibble::column_to_rownames(var = "sample"),
  annotation_colors = list(
    type =
      gdocs_pal()(7) %>%
      set_names(nm = tb_samples_extended$type %>% unique() %>% sort())
  )
)
dev.off()
