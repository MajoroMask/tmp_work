
consite <- "/data1/suna/work/tmp_work/20230825_wash_data"
# dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

R_proj <- "~/proj/tmp_work"
renv::activate(R_proj)

library(optparse)
library(vroom)
library(openxlsx)
library(rlang)
library(purrr)
library(stringr)
library(tidyr)
library(forcats)
library(dplyr)
library(fs)

# args ----

arg <-
  list(
    optparse::make_option(
      c("--pro_input"), type = "character",
      help = "Protein table input"
    ),
    optparse::make_option(
      c("--pep_input"), type = "character",
      help = "Peptide table input"
    ),
    optparse::make_option(
      c("--lc"), default = "lc_new.tsv", type = "character",
      help =
        "
        A table-like plain text file with loading control information.
        First two columns must be protein accession and mass.
        [default %default]
        "
    ),
    optparse::make_option(
      c("--output_m"), type = "character",
      help = "Output mass table"
    ),
    optparse::make_option(
      c("--output_ppm"), type = "character",
      help = "Output ppm table"
    )
  ) %>%
  optparse::OptionParser(
    usage = "Usage: %prog [options]",
    description = "Calculate protein amount based on loading control",
    option_list = .
  ) %>%
  optparse::parse_args()

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

# debug ----

arg$pro_input <- "_i/20230912/20230911-Sintilimab-1mg-Native_Digestion-STD_SC-30kd-2ug_Peptide-Direct-160min-Fix-nomod_ProtideGroups.txt"
arg$pep_input <- "_i/20230912/20230911-Sintilimab-1mg-Native_Digestion-STD_SC-30kd-2ug_Peptide-Direct-160min-Fix-nomod_PeptideGroups.txt"
arg$lc <- "_i/20230912/lc_yjj.tsv"
arg$output_m <- "_o/20230912/20230911-Sintilimab-1mg-Native_Digestion-STD_SC-30kd-2ug_Peptide-Direct-160min-Fix-nomod_m.xlsx"
arg$output_ppm <- "_o/20230912/20230911-Sintilimab-1mg-Native_Digestion-STD_SC-30kd-2ug_Peptide-Direct-160min-Fix-nomod_ppm.xlsx"

# IO ----

tb_pro <-
  path(arg$pro_input) %>%
  openxlsx::read.xlsx() %>%
  as_tibble() %>%
  select(-any_of(c("ng", "pmol")))
tb_pep <-
  path(arg$pep_input) %>%
  openxlsx::read.xlsx() %>%
  as_tibble()
# loading control
tb_lc <-
  path(arg$lc) %>%
  vroom::vroom(show_col_types = FALSE) %>%
  select(1, 2, 3) %>%
  set_names(nm = c("pro", "m", "mw"))

# main ----

tb_pro_mw <-
  tb_pro %>%
  select(
    master_pro = `Accession`,
    mw = `MW.[kDa]`
  )
tb_pep_abun <-
  tb_pep %>%
  rename(
    master_pro = `Master.Protein.Accessions`,
    anno_seq = `Annotated.Sequence`
  ) %>%
  rename_with(
    .cols = matches(stringr::regex("[Aa]bundance\\:\\.?")),
    .fn = ~ str_replace(.x, "[Aa]bundance\\:\\.?", "abun__")
  ) %>%
  filter(!str_detect(master_pro, ";")) %>%
  # filter(master_pro != "ASKG712") %>%
  select(
    master_pro,
    anno_seq,
    starts_with("abun__")
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("abun__"),
    names_to = "sample",
    values_to = "abun"
  ) %>%
  mutate(sample = str_replace(sample, "^abun__", "")) %>%
  filter(!is.na(abun)) %>%
  arrange(master_pro, sample, desc(abun))

tb_pro_tmp <-
  tb_pep_abun %>%
  filter(!master_pro %in% tb_lc$pro) %>%
  group_by(master_pro, sample) %>%
  slice_head(n = 3) %>%
  summarise(
    abun_sum = sum(abun),
    # top_n = n(),
    top_n = 3,
    .groups = "drop"
  ) %>%
  ungroup() %>%
  left_join(tb_pro_mw, by = "master_pro")

tb_lc_abun <-
  tb_pep_abun %>%
  filter(master_pro %in% tb_lc$pro) %>%
  group_by(master_pro, sample) %>%
  reframe(
    abun_sum = c(
      if_else(
        length(abun) >= 3,
        head(abun, n = 3) %>% sum(),
        NA_real_
      ),
      if_else(
        length(abun) >= 2,
        head(abun, n = 2) %>% sum(),
        NA_real_
      ),
      head(abun, n = 1)
    ),
    top_n = c(3, 2, 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(abun_sum)) %>%
  left_join(tb_lc, by = c("master_pro" = "pro")) %>%
  # left_join(tb_pro_mw, by = "master_pro") %>%
  rename(lc = master_pro) %>%
  rename_with(
    .cols = c(abun_sum, mw, m),
    .fn = ~ str_c(.x, "_lc")
  )

tb_pro_abun <-
  tb_lc_abun %>%
  left_join(
    tb_pro_tmp,
    by = c("sample", "top_n"),
    relationship = "many-to-many"
  ) %>%
  group_by(master_pro, sample) %>%
  summarise(
    m = sum(abun_sum / abun_sum_lc * m_lc / mw_lc * mw) / n(),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(sample = factor(sample)) %>%
  arrange(sample) %>%
  mutate(sample = as.character(sample)) %>%
  tidyr::pivot_wider(
    id_cols = master_pro,
    names_from = sample,
    names_prefix = "m_",
    values_from = m,
    values_fill = 0
  )
tb_pro_output <-
  tb_pro %>%
  left_join(tb_pro_abun, by = c("Accession" = "master_pro"))
openxlsx::write.xlsx(
  tb_pro_output,
  file = path(arg$output_m),
  overwrite = TRUE
)

tb_pro_output_ppm <-
  tb_pro_output %>%
  mutate(across(starts_with("m_"), .fns = ~ .x / 20)) %>%
  rename_with(
    .cols = starts_with("m_"),
    .fn = ~ str_replace(.x, "^m_", "ppm_")
  )
openxlsx::write.xlsx(
  tb_pro_output_ppm,
  file = path(arg$output_ppm),
  overwrite = TRUE
)
