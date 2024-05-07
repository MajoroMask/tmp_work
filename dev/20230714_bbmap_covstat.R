
consite <- "/data1/suna/work/tmp_work/20230714_bbmap_covstat"
dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

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

# functions ----

# magic ----

a <- new_environment()
a$pdi <- path("_i") %>% dir_create()
a$pdo <- path("_o") %>% dir_create()
# a$pd_xit <- path("/data1/wangyi/work/mg0625/") %>% path_real()
# a$pd_res <- path(a$pd_xit, "output_map_vector") %>% path_real()
a$pd_xit <- path("/data1/suna/work/mag0724/") %>% path_real()
a$pd_res <- path(a$pd_xit, "output_map_extended") %>% path_real()
a$pdo_bincov <- path(a$pdo, "bincov") %>% dir_create()
a$pdo_basecov <- path(a$pdo, "basecov") %>% dir_create()

# a$sample_order <-
#   c(
#     "293T_E5_SBHV2_E5_PPV5_E6",
#     "293T_E5_SBHV2_E5_PPV5_E6_KZ.1",
#     "293T_E5_SBHV2_E5_PPV5_E6_KZ.2",
#     "293T_E5_SBHV2_E5_PPV5_E6_QSZ_KZ.1",
#     "293T_E5_SBHV2_E5_PPV5_E6_QSZ_KZ.2",
#     "NTC.1",
#     "NTC.2"
#   )
a$sample_order <-
  c(
    "lib_236-293T_RNA_02-20230721",
    "lib_237-293T_RNA_03-20230721",
    "lib_238-293T_RNA_04-20230721",
    "lib_239-293T_RNA_05-20230721",
    "lib_240-293T_RNA_06-20230721",
    "lib_241-293T_RNA_07-20230721",
    "lib_242-293T_RNA_08-20230721",
    "lib_243-293T_RNA_09-20230721",
    "lib_244-293T_RNA_10-20230721",
    "lib_245-293T_RNA_01-20230721",
    "lib_246-NTC_01-20230721",
    "lib_247-NTC_02-20230721",
    "lib_145-293T_DNA_3D-20230427",
    "lib_173-293T_DNA_lod-20230505"
  )

# main ----

## input xit ----

tb_samples <-
  path(a$pd_xit, "samplesheet.csv") %>%
  vroom(show_col_types = FALSE) %>%
  mutate(sample = factor(sample, levels = a$sample_order)) %>%
  arrange(sample)
tb_virus <-
  path(a$pd_res, "Taxonomy", "virus_mapping", "candidate_genome_info.csv") %>%
  vroom(show_col_types = FALSE) %>%
  arrange(
    desc(in_refsheet),
    desc(in_whitelist)
  ) %>%
  mutate(tax_name = factor(tax_name) %>% fct_inorder())

## binned coverage ----

tb_files <-
  tidyr::expand_grid(
    tb_samples %>% select(sample),
    tb_virus %>% select(tax_id)
  ) %>%
  set_names(nm = c("sample", "tax_id")) %>%
  mutate(
    bincov =
      path(
        a$pd_res, "Taxonomy", "virus_mapping", sample, tax_id,
        glue::glue("bbmap_{sample}_{tax_id}_bincov.txt")
      ) %>%
      path_real(),
    basecov =
      path(
        a$pd_res, "Taxonomy", "virus_mapping", sample, tax_id,
        glue::glue("bbmap_{sample}_{tax_id}_basecov.txt")
      ) %>%
      path_real()
  )

tb_bincov <-
  purrr::pmap_dfr(
    .l = tb_files %>% select(smp = sample, tid = tax_id, file = bincov),
    .f = function(smp, tid, file) {
      tb_out <-
        vroom(
          file = file,
          skip = 3L, delim = "\t",
          col_names = c("refname", "cov", "pos", "running_pos"),
          show_col_types = FALSE
        ) %>%
        mutate(
          sample = smp,
          tax_id = tid,
          .before = 1
        ) %>%
        mutate(
          fct_refname = fct_inorder(refname) %>% as.integer(),
          .after = refname
        )
      return(tb_out)
    }
  )

tb_p_bincov <-
  tb_bincov %>%
  left_join(tb_virus %>% select(tax_id, tax_name), by = "tax_id") %>%
  # filter(tax_id == 2956715) %>%
  pipe_end()

p_bincov_overview <-
  ggplot(tb_p_bincov) +
  geom_bar(
    aes(running_pos, cov, fill = tax_name),
    stat = "identity"
  ) +
  scale_fill_tableau(palette = "Tableau 20", guide = "none") +
  facet_wrap(
    facets = vars(tax_name),
    ncol = 2, scales = "free"
  ) +
  theme_bw()

l_p_bincov <-
  purrr::map2(
    .x = tb_p_bincov %>% group_by(tax_name) %>% group_split(),
    .y = levels(tb_virus$tax_name),
    .f = function(tb, tax_name) {
      p <-
        ggplot(tb) +
        geom_bar(
          aes(running_pos, cov, fill = sample),
          stat = "identity",
          position = position_dodge()
        ) +
        scale_x_continuous(n.breaks = 10L) +
        scale_fill_tableau(
          palette = "Tableau 10",
          guide = guide_legend(
            title = NULL,
            direction = "horizontal",
            nrow = 2
          )
        ) +
        facet_wrap(
          facets = vars(tax_name),
          ncol = 1, scales = "free"
        ) +
        theme_bw() +
        theme(
          strip.background = element_rect(fill = alpha("black", 0)),
          axis.title.x = element_blank(),
          legend.position = "bottom"
        ) +
        labs(y = "Mean depth (binned per 100bp)")
      ggsave(
        filename = path(a$pdo_bincov, glue::glue("{tax_name}.pdf")),
        plot = p,
        width = 12, height = 6
      )
      return(p)
    }
  )

## 特化型basecov ----

tb_basecov <-
  purrr::pmap_dfr(
    .l =
      tb_files %>%
      select(smp = sample, tid = tax_id, file = basecov) %>%
      filter(tid %in% c(2956714, 2956715, 2956716)),
    .f = function(smp, tid, file) {
      tb_out <-
        vroom(
          file = file,
          skip = 1, delim = "\t",
          col_names = c("refname", "pos", "cov"),
          show_col_types = FALSE
        ) %>%
        mutate(
          sample = smp,
          tax_id = tid,
          .before = 1
        ) %>%
        mutate(
          fct_refname = fct_inorder(refname) %>% as.integer(),
          .after = refname
        )
      return(tb_out)
    }
  )
tb_p_basecov <-
  tb_basecov %>%
  left_join(tb_virus %>% select(tax_id, tax_name), by = "tax_id") %>%
  # filter(tax_id == 2956715) %>%
  pipe_end()

l_p_basecov <-
  purrr::map2(
    .x = tb_p_basecov %>% group_by(tax_name) %>% group_split(),
    .y = tb_p_basecov %>% group_by(tax_name) %>% group_keys() %>% pull(1),
    .f = function(tb, tax_name) {
      p <-
        tb %>%
        mutate(
          size = if_else(cov == 0, "zero", "non-zero"),
          group =
            str_c(sample, " - ", tax_name, sep = "") %>%
            factor() %>%
            fct_inorder(),
        ) %>%
        ggplot() +
        geom_point(aes(pos, cov, color = sample, size = cov)) +
        scale_x_continuous(n.breaks = 10L) +
        scale_color_tableau(
          palette = "Tableau 10",
          # guide = guide_legend(
          #   title = NULL,
          #   direction = "horizontal",
          #   nrow = 2
          # )
          guide = "none"
        ) +
        # scale_size_manual(
        #   values = c("zero" = 0.25, "non-zero" = 1.5),
        #   guide = "none"
        # ) +
        scale_size(
          range = c(0.25, 2),
          guide = "none"
        ) +
        facet_wrap(
          facets = vars(group),
          ncol = 1, scales = "free"
        ) +
        theme_bw() +
        theme(
          strip.background = element_rect(fill = alpha("black", 0)),
          axis.title.x = element_blank(),
          legend.position = "bottom"
        ) +
        labs(y = "Depth")
      ggsave(
        filename = path(a$pdo_basecov, glue::glue("basecov_{tax_name}.pdf")),
        plot = p,
        width = 8, height = 12, scale = 1.25
      )
      return(p)
    }
  )

## basecov ----

tb_basecov <-
  purrr::pmap_dfr(
    .l =
      tb_files %>%
      select(smp = sample, tid = tax_id, file = basecov) %>%
      filter(tid %in% c(11320)) %>%
      pipe_end(),
    .f = function(smp, tid, file) {
      tb_out <-
        vroom(
          file = file,
          skip = 1, delim = "\t",
          col_names = c("refname", "pos", "cov"),
          show_col_types = FALSE
        ) %>%
        mutate(
          sample = smp,
          tax_id = tid,
          .before = 1
        ) %>%
        mutate(
          fct_refname = fct_inorder(refname) %>% as.integer(),
          .after = refname
        )
      return(tb_out)
    }
  )
tb_p_basecov <-
  tb_basecov %>%
  left_join(tb_virus %>% select(tax_id, tax_name), by = "tax_id") %>%
  # filter(tax_id == 2956715) %>%
  pipe_end()

l_p_basecov <-
  purrr::map2(
    .x = tb_p_basecov %>% group_by(tax_name) %>% group_split(),
    .y = tb_p_basecov %>% group_by(tax_name) %>% group_keys() %>% pull(1),
    .f = function(tb, tax_name) {
      p <-
        tb %>%
        mutate(
          size = if_else(cov == 0, "zero", "non-zero"),
          # group =
          #   str_c(sample, " - ", tax_name, sep = "") %>%
          #   factor() %>%
          #   fct_inorder(),
          group_col = str_replace(refname, ".*(segment : \\d+$)", "\\1"),
          group_row = sample
        ) %>%
        ggplot() +
        geom_point(aes(pos, cov, color = sample, size = cov, alpha = size)) +
        scale_x_continuous(n.breaks = 10L) +
        scale_color_tableau(
          palette = "Tableau 20",
          # guide = guide_legend(
          #   title = NULL,
          #   direction = "horizontal",
          #   nrow = 2
          # )
          guide = "none"
        ) +
        # scale_size_manual(
        #   values = c("zero" = 0.25, "non-zero" = 1.5),
        #   guide = "none"
        # ) +
        scale_size(
          range = c(0.25, 2),
          guide = "none"
        ) +
        scale_alpha_manual(
          values = c("zero" = 0, "non-zero" = 1),
          guide = "none"
        ) +
        # facet_wrap(
        #   facets = vars(group),
        #   ncol = 1, scales = "free"
        # ) +
        facet_grid(
          rows = vars(group_row),
          cols = vars(group_col),
          scales = "free"
        ) +
        theme_bw() +
        theme(
          strip.background = element_rect(fill = alpha("black", 0)),
          axis.title.x = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(hjust = 1, angle = 30)
        ) +
        labs(y = "Depth")
      ggsave(
        filename = path(a$pdo_basecov, glue::glue("basecov_{tax_name}.pdf")),
        plot = p,
        width = 20, height = 15, scale = 1.25
      )
      return(p)
    }
  )
