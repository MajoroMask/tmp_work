
library(rlang)
library(vroom)
library(dplyr)
library(ggplot2)
library(fs)

setwd("~/work/tmp_work/20230427_bracken_results/")

# magic ----

the10 <- c(
  "Human immunodeficiency virus 1",
  "Bracoviriform facetosae",
  "Megavirus chiliensis",
  "Tomato spotted wilt orthotospovirus",
  "White spot syndrome virus",
  "Ichnoviriform fugitivi",
  "Caviid betaherpesvirus 2",
  "Woolly monkey sarcoma virus",
  "Choristoneura rosaceana nucleopolyhedrovirus"
)

# main ----

tb_input <-
  # path("./_i/bracken_outputs_S_ham.txt") %>%
  path("./_i/bracken_outputs_S_homo.txt") %>%
  # path("./_i/bracken_outputs_S_0504.txt") %>%
  vroom::vroom()
tb_long <-
  tb_input %>%
  tidyr::pivot_longer(
    cols = -c(name, taxonomy_id, taxonomy_lvl),
    names_pattern = "(.+)_(num|frac)",
    names_to = c("sample", "what")
  ) %>%
  mutate(
    type = stringr::str_replace(sample, "\\.\\d+", "")
    # type = case_when(
    #   sample %in% c("lib_144", "lib_147", "lib_150") ~ "CHO",
    #   sample %in% c("lib_145", "lib_148") ~ "293T",
    #   sample %in% c("lib_146", "lib_149", "lib_151", "lib_152") ~ "NTC",
    #   TRUE ~ NA_character_
    # )
  )

# view
tb_long %>%
  filter(
    type == "NTC",
    name == "Ovine mastadenovirus A"
  )
tb_long %>%
  filter(what == "num", name == "Bovine respirovirus 3") %>%
  arrange(type, sample)

enlist_name <-
  tb_long %>%
  filter(value > 50) %>%
  pull(name) %>%
  unique()

tb_p <-
  tb_long %>%
  filter(name %in% enlist_name)

p1 <-
  tb_p %>%
  filter(what == "frac") %>%
  ggplot() +
  geom_boxplot(aes(type, value, fill = type)) +
  ggsci::scale_fill_lancet() +
  facet_grid(cols = vars(what), scales = "free") +
  coord_cartesian(ylim = c(0, 0.01)) +
  labs()
p2 <-
  tb_p %>%
  filter(what == "num") %>%
  ggplot() +
  geom_boxplot(aes(type, value, fill = type)) +
  ggsci::scale_fill_lancet() +
  facet_grid(cols = vars(what), scales = "free") +
  coord_cartesian(ylim = c(0, 800)) +
  labs()
p3_num <-
  tb_p %>%
  filter(what == "num") %>%
  arrange(type) %>%
  mutate(
    x = forcats::fct_inorder(sample),
    name = forcats::fct_inorder(name)
  ) %>%
  ggplot() +
  geom_bar(
    aes(x, y = value, fill = type),
    stat = "identity",
    position = "dodge"
  ) +
  ggsci::scale_fill_lancet() +
  facet_wrap(facets = vars(name), nrow = 8, scales = "free_y") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)
  )
p3_frac <-
  tb_p %>%
  filter(what == "frac") %>%
  mutate(
    x = sample,
    name = forcats::fct_inorder(name)
  ) %>%
  ggplot() +
  geom_bar(
    aes(x, y = value, fill = type),
    stat = "identity",
    position = "dodge"
  ) +
  ggsci::scale_fill_lancet() +
  facet_wrap(facets = vars(name), nrow = 8, scales = "free_y") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)
  )
p4 <-
  tb_p %>%
  filter(
    what == "num",
    # type == "293T_BAV2"
    type == "293T"
  ) %>%
  arrange(desc(value)) %>%
  mutate(name = forcats::fct_inorder(name)) %>%
  ggplot() +
  geom_bar(
    aes(name, value + 1, fill = sample),
    stat = "identity",
    position = "dodge"
  ) +
  scale_y_log10() +
  ggsci::scale_fill_lancet() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)
  )
p5 <-
  tb_p %>%
  filter(
    what == "num",
    # type == "BAV2"
    type == "CHO"
  ) %>%
  arrange(desc(value)) %>%
  mutate(name = forcats::fct_inorder(name)) %>%
  ggplot() +
  geom_bar(
    aes(name, value + 1, fill = sample),
    stat = "identity",
    position = "dodge"
  ) +
  scale_y_log10() +
  ggsci::scale_fill_lancet() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)
  )
p6 <-
  tb_p %>%
  filter(
    what == "num",
    type == "NTC"
  ) %>%
  arrange(desc(value)) %>%
  mutate(name = forcats::fct_inorder(name)) %>%
  ggplot() +
  geom_bar(
    aes(name, value + 1, fill = sample),
    stat = "identity",
    position = "dodge"
  ) +
  scale_y_log10() +
  ggsci::scale_fill_lancet() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

# update 2023-05-04 ----

# 不做host removal的话，293T数据中的Bovine respirovirus 3是这样的
# btw，数据直接从原始数据文件里copy进来

sample_names <-
  fs::path("./_i/bracken") %>%
  fs::dir_ls() %>%
  fs::path_file()
tb_nohr <-
  purrr::map(
    sample_names,
    .f = ~
      glue::glue("./_i/bracken/{.x}/{.x}_S.bracken_output.txt") %>%
      fs::path() %>%
      vroom::vroom(show_col_types = FALSE) %>%
      select(
        name, taxonomy_id, taxonomy_lvl,
        new_est_reads, fraction_total_reads
      ) %>%
      rename(
        "{.x}_num" := new_est_reads,
        "{.x}_frac" := fraction_total_reads,
      )
  ) %>%
  purrr::reduce(
    .f = full_join,
    by = c("name", "taxonomy_id", "taxonomy_lvl")
  ) %>%
  mutate(
    across(
      .col = -c(name, taxonomy_id, taxonomy_lvl),
      .fns = ~ tidyr::replace_na(.x, replace = 0)
    )
  )
tb_hr <-
  path("./_i/bracken_outputs_S_homo.txt") %>%
  vroom::vroom(show_col_types = FALSE) %>%
  select(
    name, taxonomy_id, taxonomy_lvl,
    all_of(glue::glue("{sample_names}_num"))
  ) %>%
  filter(`293T_BAV2.1_num` > 0 | `293T_BAV2.2_num` > 0 | `293T_BAV2.3_num` > 0)

# 一些统计

VennDiagram::venn.diagram(
  x = list(
    hr = tb_hr$name,
    no_hr = tb_nohr$name
  ),
  filename = "./_o/venn.pdf"
)

tmp_tb2m <- function(tb) {
  # s_in_both <- intersect(tb_hr$name, tb_nohr$name)
  sample_names <- c("293T_BAV2.1", "293T_BAV2.2", "293T_BAV2.3")
  those_columns <- c("name", glue::glue("{sample_names}_num"))

  m <-
    tb %>%
    # filter(name %in% s_in_both) %>%
    select(all_of(those_columns)) %>%
    arrange(name) %>%
    tibble::column_to_rownames(var = "name") %>%
    as.matrix()
  return(m)
}

# tmp_cosine_simi <- function(x, y) {
#   x %*% y / sqrt(x%*%x * y%*%y)
# }

s_in_both <- intersect(tb_hr$name, tb_nohr$name)

m_hr_all <- tmp_tb2m(tb_hr)
m_nohr_all <- tmp_tb2m(tb_nohr)

# 相似性和倍数

m_hr <- m_hr_all[rownames(m_hr_all) %in% s_in_both, ]
m_nohr <- m_nohr_all[rownames(m_nohr_all) %in% s_in_both, ]
cor(m_hr, m_nohr)
lsa::cosine(cbind(m_hr, m_nohr))
norm(m_hr, type = "F") / norm(m_nohr, type = "F")
sum(m_hr) / sum(m_nohr)

#

m_hr_uniq <- m_hr_all[!rownames(m_hr_all) %in% s_in_both, ]
m_nohr_uniq <- m_nohr_all[!rownames(m_nohr_all) %in% s_in_both, ]
summary(m_hr_uniq)
summary(m_nohr_uniq)
summary(m_nohr)

tb_comb <-
  full_join(
    tb_hr, tb_nohr,
    by = c("name", "taxonomy_id", "taxonomy_lvl"),
    suffix = c("_hr", "_nohr")
  ) %>%
  tidyr::pivot_longer(
    cols = -c(name, taxonomy_id, taxonomy_lvl),
    names_pattern = "(.+)_(num|frac)_(hr|nohr)",
    names_to = c("sample", "what", "source")
  ) %>%
  filter(sample %in% sample_names)

# xit <-
#   tribble(
#     ~name,                   ~tax_id, ~tax_lvl, ~sample,       ~kraken_assigned_reads, ~added_reads, ~new_est_reads, ~fraction_total_reads, ~source,
#     "Bovine respirovirus 3", 11215,   "S",      "293T_BAV2.1", 4,                      0,            4,              0.00000,               "no host removal",
#     "Bovine respirovirus 3", 11215,   "S",      "293T_BAV2.2", 2,                      0,            2,              0.00000,               "no host removal",
#     "Bovine respirovirus 3", 11215,   "S",      "293T_BAV2.3", 8,                      0,            8,              0.00000,               "no host removal",
#     "Bovine respirovirus 3", 11215,   "S",      "293T_BAV2.1", 1,                      0,            1,              0.00000,               "host removal",
#     "Bovine respirovirus 3", 11215,   "S",      "293T_BAV2.2", 0,                      0,            0,              0.00000,               "host removal",
#     "Bovine respirovirus 3", 11215,   "S",      "293T_BAV2.3", 1,                      0,            1,              0.00000,               "host removal"
#   ) %>%
#   select(-tax_lvl)

# update 2023-05-05 ----

# 给明明合并下数据

tb_input <-
  path("./_i/bracken_outputs_S_0504.txt") %>%
  vroom::vroom(show_col_types = FALSE)
tb_from_hmm <-
  path("./_i/3D_Virus_counts.xlsx") %>%
  openxlsx::read.xlsx(sheet = 1) %>%
  tibble::as_tibble() %>%
  slice(-1)
tb_comb <-
  full_join(
    tb_from_hmm,
    tb_input %>%
      select(name, ends_with("_num")) %>%
      rename_with(
        .cols = -name,
        .fn = ~ stringr::str_replace(.x, "_num$", "")
      ),
    by = c("Name" = "name")
  ) %>%
  mutate(
    across(
      .cols = -Name,
      .fns = ~ tidyr::replace_na(.x, 0)
    )
  )
openxlsx::write.xlsx(
  tb_comb,
  file = "./_o/combined_sheet.xlsx",
  overwrite = TRUE
)
