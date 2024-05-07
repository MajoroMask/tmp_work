
consite <- "/data1/suna/work/tmp_work/20231206_plot_venn"
dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

library(openxlsx)
library(ggVennDiagram)

library(rlang)
library(glue)
library(purrr)
library(stringr)
library(tidyr)
library(dplyr)
library(fs)

# magic ----

a <- new_environment()
a$pdi <- path("_i") %>% dir_create()
a$pdo <- path("_o") %>% dir_create()

# func ----

get_blennded_colors <- function(colors, try_alpha = c(0.5, 1)) {
  n_ele <- length(colors)
  alpha_step <- seq(try_alpha[1], try_alpha[2], by = 1/n_ele)[-1]

  tb_xit <-
    rep(list(c(0, 1)), n_ele) %>%
    set_names(nm = paste0("xit", 1:n_ele)) %>%
    as.data.frame() %>%
    expand.grid() %>%
    as_tibble() %>%
    slice(-1) %>%
    rowwise() %>%
    mutate(
      id =
        paste0(
          colnames(.)[as.logical(c_across(cols = everything()))],
          collapse = ""
        ) %>%
        str_replace_all("xit", ""),
      .sum = sum(c_across(-id)),
      blend_alpha = alpha_step[.sum],
      .before = 1
    ) %>%
    ungroup() %>%
    arrange(.sum, id) %>%
    select(-.sum) %>%
    rowwise() %>%
    mutate(
      blend_color = colorjam::blend_colors(
        x = colors[as.logical(c_across(cols = starts_with("xit")))]
      )
    ) %>%
    select(id, blend_color, blend_alpha)
  return(tb_xit)
}

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

# magic ----

the_five <- c(
  "GOI4-MS-P1-F-Par",
  "GOI4-MS-P2-F-Par",
  "GOI4-MS-P3-F-Par",
  "GOI4-MS-S-F-Par",
  "GOI4-MS-E-F-Par"
)

# main ----

## venn diagram ----

### data.xlsx ----

{
  wb <- loadWorkbook(xlsxFile = path(a$pdi, "data.xlsx"))
  preset_colnames <-
    c("acc", "desc", "mark", "n_unipep", "n_aa", "mw", "pi", "m_mean")
  l_data_full <- purrr::map(
    .x = names(wb),
    .f = ~
      readWorkbook(wb, sheet = .x) %>%
      set_names(nm = preset_colnames)
  )
  names(l_data_full) <- names(wb)

  l_data_sub <- purrr::map(
    .x = l_data_full,
    .f = ~
      .x %>%
      mutate(pct_m = m_mean / sum(m_mean)) %>%
      filter(pct_m >= 0.001)
  )
  names(l_data_sub) <- names(wb)
}

### data_full.xlsx ----

{
  wb <- loadWorkbook(xlsxFile = path(a$pdi, "data_p2.xlsx"))
  preset_colnames <-
    c("acc", "desc", "n_unipep", "n_aa", "mw", "pi", "m_mean")
  l_data_p2 <- purrr::map(
    .x = names(wb)[-1],
    .f = ~
      readWorkbook(wb, sheet = .x) %>%
      set_names(nm = preset_colnames) %>%
      select(acc, m_mean)
  )
  names(l_data_p2) <- names(wb)[-1]

  l_data_p2_sub <- purrr::map(
    .x = l_data_p2,
    .f = ~
      .x %>%
      mutate(pct_m = m_mean / sum(m_mean)) %>%
      filter(pct_m >= 0.001)
  )
  names(l_data_p2_sub) <- names(wb)[-1]
}

### plot pentagon ----

for (i in 0:8) {
  l_xit <- l_data_sub
  if (i != 0) {
    l_xit$`GOI4-MS-P2-F-Par` <- l_data_p2_sub[[i]]
  }
  l_data <- l_xit %>% map(.f = ~ .x %>% pull(acc))

  venn_penta <- Venn(l_data[the_five])
  data_penta <- process_data(venn_penta)

  p_penta <-
    ggplot(
      data =
        venn_region(data_penta) %>%
        left_join(
          get_blennded_colors(colors = ggthemes::gdocs_pal()(5)),
          by = "id"
        )
    ) +
    geom_sf(aes(fill = blend_color, alpha = blend_alpha)) +
    geom_sf(
      size = 2,
      linetype = "dashed",
      color = "grey",
      data = venn_setedge(data_penta),
      show.legend = TRUE
    ) +
    geom_sf_text(aes(label = name), data = venn_setlabel(data_penta)) +
    geom_sf_label(
      aes(
        label = ifelse(
          count == 0,
          count,
          paste0(
            count, " (",
            scales::percent(count / sum(count), accuracy = 0.1),
            ")"
          )
        )
      )
    ) +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    scale_fill_identity() +
    scale_alpha_identity() +
    theme_void()

  ggsave(
    p_penta, filename = path(a$pdo, glue("penta_{i}.pdf")),
    width = 10, height = 8
  )
  ggsave(
    p_penta, filename = path(a$pdo, glue("penta_{i}.tiff")),
    width = 5, height = 4
  )
}

for (i in 0:8) {
  l_xit <- l_data_full
  if (i != 0) {
    l_xit$`GOI4-MS-P2-F-Par` <- l_data_p2[[i]]
  }
  l_data <- l_xit %>% map(.f = ~ .x %>% pull(acc))

  venn_penta <- Venn(l_data[the_five])
  data_penta <- process_data(venn_penta)

  p_penta <-
    ggplot(
      data =
        venn_region(data_penta) %>%
        left_join(
          get_blennded_colors(colors = ggthemes::gdocs_pal()(5)),
          by = "id"
        )
    ) +
    geom_sf(aes(fill = blend_color, alpha = blend_alpha)) +
    geom_sf(
      size = 2,
      linetype = "dashed",
      color = "grey",
      data = venn_setedge(data_penta),
      show.legend = TRUE
    ) +
    geom_sf_text(aes(label = name), data = venn_setlabel(data_penta)) +
    geom_sf_label(
      aes(
        label = ifelse(
          count == 0,
          count,
          paste0(
            count, " (",
            scales::percent(count / sum(count), accuracy = 0.1),
            ")"
          )
        )
      )
    ) +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    scale_fill_identity() +
    scale_alpha_identity() +
    theme_void()

  ggsave(
    p_penta, filename = path(a$pdo, glue("penta_all_{i}.pdf")),
    width = 10, height = 8
  )
  ggsave(
    p_penta, filename = path(a$pdo, glue("penta_all_{i}.tiff")),
    width = 5, height = 4
  )
}

### plot tri ----

l_the_three <- list(
  c("GOI4-MS-P1-F-Par", "GOI4-MS-P1-F-Sup", "GOI4-MS-P1-F-Mix"),
  c("GOI4-MS-P3-F-Par", "GOI4-MS-P3-F-Sup", "GOI4-MS-P3-F-Mix")
)

for (i in 1:2) {
  the_three <- l_the_three[[i]]
  the_name <- c("P1", "P3")[i]
  l_data <- l_data_full %>% map(.f = ~ .x %>% pull(acc))

  venn_tri <- Venn(l_data[the_three])
  data_tri <- process_data(venn_tri)

  p_penta <-
    ggplot(
      data =
        venn_region(data_tri) %>%
        left_join(
          get_blennded_colors(colors = ggthemes::gdocs_pal()(3)),
          by = "id"
        )
    ) +
    geom_sf(aes(fill = blend_color, alpha = blend_alpha)) +
    geom_sf(
      size = 2,
      linetype = "dashed",
      color = "grey",
      data = venn_setedge(data_tri),
      show.legend = TRUE
    ) +
    geom_sf_text(aes(label = name), data = venn_setlabel(data_tri)) +
    geom_sf_label(
      aes(
        label = ifelse(
          count == 0,
          count,
          paste0(
            count, " (",
            scales::percent(count / sum(count), accuracy = 0.1),
            ")"
          )
        )
      )
    ) +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    scale_fill_identity() +
    scale_alpha_identity() +
    theme_void()

  ggsave(
    p_penta, filename = path(a$pdo, glue("tri_{the_name}.pdf")),
    width = 10, height = 8
  )
  ggsave(
    p_penta, filename = path(a$pdo, glue("tri_{the_name}.tiff")),
    width = 5, height = 4
  )
}

### cor ----

the_three <- c("GOI4-MS-P1-F-Par", "GOI4-MS-P1-F-Sup", "GOI4-MS-P1-F-Mix")
# the_three <- c("GOI4-MS-P3-F-Par", "GOI4-MS-P3-F-Sup", "GOI4-MS-P3-F-Mix")

tb_cov <-
  tb_data_ori %>%
  filter(sample %in% the_three) %>%
  mutate(
    sample_short = str_replace(sample, ".*-", ""),
    .after = sample
  ) %>%
  # select(sample_short, acc, m_mean) %>%
  select(sample_short, acc, pct_m, in_top_pct) %>%
  pivot_wider(
    id_cols = acc,
    names_from = sample_short,
    # values_from = m_mean
    values_from = c(pct_m, in_top_pct)
  ) %>%
  rowwise() %>%
  mutate(any_sample = any(c_across(starts_with("in_top_pct_")))) %>%
  ungroup() %>%
  filter(any_sample) %>%
  pipe_end()
m_cov <-
  tb_cov %>%
  select(acc, starts_with("pct_m_")) %>%
  rename_with(
    .cols = starts_with("pct_m_"),
    .fn = ~ str_replace(.x, "pct_m_", "")
  ) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "acc") %>%
  as.matrix() %>%
  na.omit()
cor(m_cov)


## group analysis ----

### data ----

{
  wb <- loadWorkbook(xlsxFile = path(a$pdi, "data.xlsx"))
  preset_colnames <-
    c("acc", "desc", "mark", "n_unipep", "n_aa", "mw", "pi", "m_mean")
  tb_data_ori <-
    purrr::map_dfr(
      .x = names(wb),
      .f = ~
        readWorkbook(wb, sheet = .x) %>%
        set_names(nm = preset_colnames) %>%
        as_tibble() %>%
        mutate(sample = .x, .before = 1) %>%
        mutate(
          pct_m = m_mean / sum(m_mean),
          in_top_pct = pct_m >= 0.01
          # in_top_pct = pct_m >= 0.01
        )
    )
}

## cor ----

the_three <- c("GOI4-MS-P1-F-Par", "GOI4-MS-P1-F-Sup", "GOI4-MS-P1-F-Mix")
# the_three <- c("GOI4-MS-P3-F-Par", "GOI4-MS-P3-F-Sup", "GOI4-MS-P3-F-Mix")

tb_cov <-
  tb_data_ori %>%
  filter(sample %in% the_three) %>%
  mutate(
    sample_short = str_replace(sample, ".*-", ""),
    .after = sample
  ) %>%
  # select(sample_short, acc, m_mean) %>%
  select(sample_short, acc, pct_m, in_top_pct) %>%
  pivot_wider(
    id_cols = acc,
    names_from = sample_short,
    # values_from = m_mean
    values_from = c(pct_m, in_top_pct)
  ) %>%
  rowwise() %>%
  mutate(any_sample = any(c_across(starts_with("in_top_pct_")))) %>%
  ungroup() %>%
  filter(any_sample) %>%
  pipe_end()
m_cov <-
  tb_cov %>%
  select(acc, starts_with("pct_m_")) %>%
  rename_with(
    .cols = starts_with("pct_m_"),
    .fn = ~ str_replace(.x, "pct_m_", "")
  ) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "acc") %>%
  as.matrix() %>%
  na.omit()
cor(m_cov)

### update: cor on penta ----

tb_cov <-
  tb_data_ori %>%
  filter(sample %in% the_five) %>%
  mutate(
    sample_short = str_replace(sample, "GOI4-MS-(.*)-F-Par", "\\1"),
    .after = sample
  ) %>%
  mutate(in_top_pct = pct_m >= 0) %>%
  mutate(value = pct_m) %>%
  select(sample_short, acc, value, in_top_pct) %>%
  pivot_wider(
    id_cols = acc,
    names_from = sample_short,
    values_from = c(value, in_top_pct),
    values_fill = 0
  ) %>%
  rowwise() %>%
  mutate(
    any_sample =
      any(c_across(starts_with("in_top_pct_")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(any_sample) %>%
  pipe_end()
m_cov <-
  tb_cov %>%
  select(acc, starts_with("value_")) %>%
  rename_with(
    .cols = starts_with("value_"),
    .fn = ~ str_replace(.x, "value_", "")
  ) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "acc") %>%
  as.matrix() %>%
  na.omit()
cor(m_cov)

readr::write_tsv(
  cor(m_cov) %>% as.data.frame() %>% tibble::rownames_to_column(var = "acc"),
  file = path(a$pdo, "m_cov.tsv")
)

## heatmap ----

the_three <- c("GOI4-MS-P1-F-Par", "GOI4-MS-P1-F-Sup", "GOI4-MS-P1-F-Mix")

tb_heatmap <-
  tb_data_ori %>%
  filter(sample %in% the_three) %>%
  mutate(
    sample_short = str_replace(sample, ".*-", ""),
    .after = sample
  ) %>%
  rename(value = pct_m) %>%
  # mutate(value = log10(value)) %>%
  select(sample_short, acc, value, in_top_pct) %>%
  pivot_wider(
    id_cols = acc,
    names_from = sample_short,
    values_from = c(value, in_top_pct)
  ) %>%
  rowwise() %>%
  mutate(any_sample = any(c_across(starts_with("in_top_pct_")))) %>%
  ungroup() %>%
  filter(any_sample) %>%
  pipe_end()
m_heatmap <-
  tb_heatmap %>%
  select(acc, starts_with("value_")) %>%
  rename_with(
    .cols = starts_with("value_"),
    .fn = ~ str_replace(.x, "value_", "")

  ) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "acc") %>%
  as.matrix() %>%
  na.omit() %>%
  pipe_end()
pheatmap::pheatmap(
  m_heatmap
)

## barplot first take ----

# the_three <- c("GOI4-MS-P1-F-Par", "GOI4-MS-P1-F-Sup", "GOI4-MS-P1-F-Mix")
the_three <- c("GOI4-MS-P3-F-Par", "GOI4-MS-P3-F-Sup", "GOI4-MS-P3-F-Mix")
the_seven_pros <-
  c("Gag", "Integrase", "Protease", "Rev", "Reverse", "RNase", "VSV-G")

tb_barplot <-
  tb_data_ori %>%
  filter(
    sample %in% the_three,
    acc %in% the_seven_pros
  ) %>%
  mutate(
    sample_short = str_replace(sample, ".*-", ""),
    .after = sample
  ) %>%
  mutate(value = as.integer(pct_m * 100000)) %>%
  select(sample_short, acc, value, in_top_pct) %>%
  # pivot_wider(
  #   id_cols = acc,
  #   names_from = sample_short,
  #   values_from = c(value, in_top_pct)
  # ) %>%
  # rowwise() %>%
  # mutate(any_sample = any(c_across(starts_with("in_top_pct_")))) %>%
  # ungroup() %>%
  # filter(any_sample) %>%
  pipe_end()

p_bar <- ggstatsplot::ggbarstats(
  tb_barplot,
  x = acc, y = sample_short,
  count = value,
  label = "both",
  package = "ggthemes",
  palette = "Tableau_10"
)
# ggsave(p_bar, filename = path(a$pdo, "bar_P1.pdf"), height = 10, width = 8)
ggsave(p_bar, filename = path(a$pdo, "bar_P3.pdf"), height = 10, width = 8)

## barplot ----

# the_three <- c("GOI4-MS-P1-F-Par", "GOI4-MS-P1-F-Sup", "GOI4-MS-P1-F-Mix")
the_three <- c("GOI4-MS-P3-F-Par", "GOI4-MS-P3-F-Sup", "GOI4-MS-P3-F-Mix")

the_seven_pros <-
  c("Gag", "Integrase", "Protease", "Rev", "Reverse", "RNase", "VSV-G")
enlist_pros <-
  tb_data_ori %>%
  filter(
    sample %in% the_three,
    pct_m >= 0.01
  ) %>%
  pull(acc) %>%
  unique() %>%
  # union(the_seven_pros) %>%
  pipe_end()

tb_barplot <-
  tb_data_ori %>%
  filter(sample %in% the_three) %>%
  group_by(sample) %>%
  group_map(
    .f = function(tb, id) {
      tb_out <-
        tb %>%
        mutate(tmp_grp = if_else(acc %in% enlist_pros, acc, "other")) %>%
        group_by(tmp_grp) %>%
        summarise(
          pct_m = sum(pct_m),
          sample = sample[1]
        ) %>%
        ungroup() %>%
        rename(acc = tmp_grp) %>%
        pipe_end()
    },
    .keep = TRUE
  ) %>%
  bind_rows() %>%
  mutate(
    sample_short = str_replace(sample, ".*-", ""),
    .after = sample
  ) %>%
  mutate(value = as.integer(pct_m * 100000)) %>%
  arrange(desc(value)) %>%
  mutate(acc = fct_inorder(acc)) %>%
  select(sample_short, acc, value) %>%
  pipe_end()

p_bar <-
  ggstatsplot::ggbarstats(
    tb_barplot,
    x = acc, y = sample_short,
    counts = value,
    type = "bayes",
    results.subtitle = FALSE,
    perc.k = 1,
    label = "percentage",
    package = "ggthemes",
    palette = "Tableau_20",
    legend.title = "Protein",
    xlab = "P3"
  )
# p_bar
# ggsave(p_bar, filename = path(a$pdo, "bar_P1.pdf"), height = 10, width = 8)
ggsave(p_bar, filename = path(a$pdo, "bar_P3.pdf"), height = 10, width = 8)

## similarity between P1 & P3 ----

p1_samples <- c("GOI4-MS-P1-F-Par", "GOI4-MS-P1-F-Sup", "GOI4-MS-P1-F-Mix")
p3_samples <- c("GOI4-MS-P3-F-Par", "GOI4-MS-P3-F-Sup", "GOI4-MS-P3-F-Mix")
enlist_pros <-
  tb_data_ori %>%
  filter(
    sample %in% union(p1_samples, p3_samples),
    pct_m >= 0.01
  ) %>%
  pull(acc) %>%
  unique() %>%
  # union(the_seven_pros) %>%
  pipe_end()

m_p1 <-
  tb_data_ori %>%
  filter(
    sample %in% p1_samples,
    acc %in% enlist_pros
  ) %>%
  mutate(
    sample_short = str_replace(sample, ".*-", ""),
    .after = sample
  ) %>%
  rename(value = pct_m) %>%
  select(sample_short, acc, value) %>%
  pivot_wider(
    id_cols = acc,
    names_from = sample_short,
    values_from = value
  ) %>%
  arrange(acc) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "acc") %>%
  as.matrix() %>%
  pipe_end()
m_p3 <-
  tb_data_ori %>%
  filter(
    sample %in% p3_samples,
    acc %in% enlist_pros
  ) %>%
  mutate(
    sample_short = str_replace(sample, ".*-", ""),
    .after = sample
  ) %>%
  rename(value = pct_m) %>%
  select(sample_short, acc, value) %>%
  pivot_wider(
    id_cols = acc,
    names_from = sample_short,
    values_from = value
  ) %>%
  arrange(acc) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "acc") %>%
  as.matrix() %>%
  pipe_end()
cor(as.vector(m_p1), as.vector(m_p3))

## facetted bar ----

the_seven_pros <-
  c("Gag", "Integrase", "Protease", "Rev", "Reverse", "RNase", "VSV-G")
p1_samples <- c("GOI4-MS-P1-F-Par", "GOI4-MS-P1-F-Sup", "GOI4-MS-P1-F-Mix")
p3_samples <- c("GOI4-MS-P3-F-Par", "GOI4-MS-P3-F-Sup", "GOI4-MS-P3-F-Mix")
enlist_pros <-
  tb_data_ori %>%
  filter(
    sample %in% union(p1_samples, p3_samples),
    pct_m >= 0.01
  ) %>%
  pull(acc) %>%
  unique() %>%
  union(the_seven_pros) %>%
  pipe_end()

tb_facet_bar <-
  tb_data_ori %>%
  filter(
    sample %in% union(p1_samples, p3_samples),
    acc %in% enlist_pros
  ) %>%
  mutate(
    sample_short = str_replace(sample, ".*-", ""),
    grp = if_else(sample %in% p1_samples, "P1", "P3")
  ) %>%
  arrange(desc(pct_m)) %>%
  mutate(acc = fct_inorder(acc))

p_facet_bar_the_seven <-
  ggplot(tb_facet_bar %>% filter(acc %in% the_seven_pros)) +
  geom_bar(
    aes(x = pct_m, y = sample_short, fill = sample_short),
    stat = "identity"
  ) +
  scale_x_continuous(labels = scales::percent) +
  ggthemes::scale_fill_gdocs(guide = "none") +
  facet_grid(
    rows = vars(grp),
    cols = vars(acc),
    scales = "free"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = NULL, y = NULL)
ggsave(
  p_facet_bar_the_seven,
  filename = path(a$pdo, "facet_bar_the_seven.pdf"),
  width = 18, height = 4
)

p_facet_bar_others <-
  ggplot(tb_facet_bar %>% filter(!acc %in% the_seven_pros)) +
  geom_bar(
    aes(x = pct_m, y = sample_short, fill = sample_short),
    stat = "identity"
  ) +
  scale_x_continuous(labels = scales::percent) +
  ggthemes::scale_fill_gdocs(guide = "none") +
  facet_grid(
    rows = vars(grp),
    cols = vars(acc),
    scales = "free"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = NULL, y = NULL)
ggsave(
  p_facet_bar_others,
  filename = path(a$pdo, "facet_bar_others.pdf"),
  width = 20, height = 4
)

### update: facet bar for the five ----

the_seven_pros <-
  c("Gag", "Integrase", "Protease", "Rev", "Reverse", "RNase", "VSV-G")
enlist_pros <-
  tb_data_ori %>%
  filter(
    sample %in% the_five,
    pct_m >= 0.01
  ) %>%
  pull(acc) %>%
  unique() %>%
  union(the_seven_pros) %>%
  pipe_end()

tb_facet_bar <-
  tb_data_ori %>%
  filter(
    sample %in% the_five,
    acc %in% enlist_pros
  ) %>%
  mutate(
    sample_short =
      str_replace(sample, "GOI4-MS-(.*)-F-Par", "\\1") %>%
      fct_rev(),
    grp = sample_short
  ) %>%
  arrange(desc(pct_m)) %>%
  mutate(acc = fct_inorder(acc))

p_facet_bar_the_seven <-
  ggplot(tb_facet_bar %>% filter(acc %in% the_seven_pros)) +
  geom_bar(
    aes(x = pct_m, y = sample_short, fill = sample_short),
    stat = "identity"
  ) +
  scale_x_continuous(labels = scales::percent) +
  ggthemes::scale_fill_gdocs(guide = "none") +
  facet_grid(
    # rows = vars(grp),
    cols = vars(acc),
    scales = "free"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = NULL, y = NULL)
ggsave(
  p_facet_bar_the_seven,
  filename = path(a$pdo, "facet_bar_the_seven_penta.pdf"),
  width = 18, height = 4
)

p_facet_bar_others <-
  ggplot(tb_facet_bar %>% filter(!acc %in% the_seven_pros)) +
  geom_bar(
    aes(x = pct_m, y = sample_short, fill = sample_short),
    stat = "identity"
  ) +
  scale_x_continuous(labels = scales::percent) +
  ggthemes::scale_fill_gdocs(guide = "none") +
  facet_grid(
    # rows = vars(grp),
    cols = vars(acc),
    scales = "free"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = NULL, y = NULL)
ggsave(
  p_facet_bar_others,
  filename = path(a$pdo, "facet_bar_others_penta.pdf"),
  width = 20, height = 4
)
