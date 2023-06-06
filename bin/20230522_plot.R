
# init ----

library(ggstatsplot)
library(openxlsx)
library(rlang)
library(readr)
library(vroom)
library(stringr)
library(dplyr)
# library(ggplot2)
library(fs)

# library(UniProt.ws)

pwd <- "/data1/suna/work/tmp_work/20230522_plot" %>% dir_create()
setwd("/data1/suna/work/tmp_work/20230522_plot")

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

# 1 download data ----

tb_anno <- vroom::vroom(file = "xits.tsv")

# cc_pool <-
#   tb_anno$`Subcellular location [CC]` %>%
#   na.omit() %>%
#   str_replace("^SUBCELLULAR LOCATION: ", "") %>%
#   str_replace("\\.$", "") %>%
#   str_replace("Note=.*", "") %>%
#   str_split(pattern = "\\.\\s?|\\;\\s?") %>%
#   unlist() %>%
#   unique() %>%
#   sort()
# cc_pool_only_desc <-
#   cc_pool %>%
#   str_replace(" \\{.*", "") %>%
#   unique() %>%
#   sort()
# writeLines(cc_pool, con = "cc_pool.txt")
# writeLines(cc_pool_only_desc, con = "cc_pool_only_desc.txt")

tb_results <-
  tb_anno %>%
  transmute(
    entry = Entry,
    protein_names = `Protein names`,
    gene_name = `Gene Names`,
    go_cc = `Gene Ontology (cellular component)`,
    subcellular = `Subcellular location [CC]`,
    transmembrane = Transmembrane
  ) %>%
  mutate(
    cc_has_membrane =
      str_detect(go_cc, regex("membrane", ignore_case = TRUE)) %>%
      tidyr::replace_na(replace = FALSE),
    sl_has_membrane =
      str_detect(subcellular, regex("membrane", ignore_case = TRUE)) %>%
      tidyr::replace_na(replace = FALSE),
    has_trans = !is.na(transmembrane),
    combine1 = cc_has_membrane | sl_has_membrane | has_trans,
    combine2 = combine1 & (sl_has_membrane | has_trans),
    combine3 = cc_has_membrane | sl_has_membrane,
    combine4 = has_trans
  )

# x <-
#   tb_results %>%
#   filter(cc_contain_membrane & !sl_contain_membrane & !has_trans) %>%
#   select(entry) %>%
#   left_join(tb_anno, by = c("entry" = "Entry"))
# View(x)

# x <-
#   tb_anno %>%
#   filter(str_detect(`Subcellular location [CC]`, fixed("{ECO:0000250|UniProtKB:P32455}")))
# x

# 2 plot from PD output ----

tb_input <-
  openxlsx::read.xlsx("20230522_input.xlsx") %>%
  as_tibble()
tb_p <-
  tb_input %>%
  left_join(tb_results, by = c("Accession" = "entry")) %>%
  rename(
    abundance1 = `Abundance.Ratio:.(vero-ccf)./.(virus)`,
    abundance2 = `Abundance.Ratio:.(vero-cell)./.(virus)`,
    MW = `MW.[kDa]`,
    pI = `calc..pI`
  ) %>%
  mutate(
    log10_MW = log10(MW),
    log10_abun1 = log10(abundance1) %>% na_if(2) %>% na_if(-2),
    log10_abun2 = log10(abundance2) %>% na_if(2) %>% na_if(-2),
    log2_abun1 = na_if(abundance1, 100) %>% na_if(0.01) %>% log2() ,
    log2_abun2 = na_if(abundance2, 100) %>% na_if(0.01) %>% log2() ,
    is_membrane1 =
      case_when(combine1 ~ "yes", !combine1 ~ "no") %>%
      factor(levels = c("yes", "no")),
    is_membrane2 =
      case_when(combine2 ~ "yes", !combine2 ~ "no") %>%
      factor(levels = c("yes", "no")),
    is_membrane3 =
      case_when(combine3 ~ "yes", !combine3 ~ "no") %>%
      factor(levels = c("yes", "no")),
    is_membrane4 =
      case_when(combine4 ~ "yes", !combine4 ~ "no") %>%
      factor(levels = c("yes", "no")),
  )

p_mw <- ggstatsplot::ggbetweenstats(
  tb_p,
  x = is_membrane1,
  y = log10_MW,
  xlab = "Is membrane protein",
  ylab = "Molecular weight (log10 scaled)",
  bf.message = FALSE
)
p_pi <- ggstatsplot::ggbetweenstats(
  tb_p,
  x = is_membrane1,
  y = pI,
  xlab = "Is membrane protein",
  ylab = "pI",
  bf.message = FALSE
)
p_ab1 <- ggstatsplot::ggbetweenstats(
  tb_p,
  x = is_membrane1,
  y = log2_abun1,
  xlab = "Is membrane protein",
  ylab = "Abundance vero-cff/virus (log2 scaled)",
  bf.message = FALSE,
)
p_ab2 <- ggstatsplot::ggbetweenstats(
  tb_p,
  x = is_membrane1,
  y = log2_abun2,
  xlab = "Is membrane protein",
  ylab = "Abundance vero-cell/virus (log2 scaled)",
  bf.message = FALSE,
)

# 3 output ----

tb_output <-
  tb_input %>%
  left_join(tb_results, by = c("Accession" = "entry")) %>%
  mutate(
    ccf_in_range =
      `Abundance.Ratio:.(vero-ccf)./.(virus)` >= 0.25 &
      `Abundance.Ratio:.(vero-ccf)./.(virus)` <= 4,
    cell_in_range =
      `Abundance.Ratio:.(vero-cell)./.(virus)` >= 0.25 &
      `Abundance.Ratio:.(vero-cell)./.(virus)` <= 4
  )
openxlsx::write.xlsx(tb_output, "20230522_output.xlsx", overwrite = TRUE)

tb_output %>%
  gsa(ccf_in_range, combine1) %>%
  readr::write_tsv("ccf_combine1.tsv")
tb_output %>%
  gsa(cell_in_range, combine1) %>%
  readr::write_tsv("cell_combine1.tsv")

ggsave(p_mw, file = "p_mw.pdf", height = 8, width = 6)
ggsave(p_pi, file = "p_pi.pdf", height = 8, width = 6)
ggsave(p_ab1, file = "p_ab1.pdf", height = 8, width = 6)
ggsave(p_ab2, file = "p_ab2.pdf", height = 8, width = 6)

# 4 another plot ----

tb_293T_in <- vroom::vroom("20230519-DDA-293T-120ug-protein.csv")
tb_p2 <-
  tb_293T_in %>%
  mutate(
    pi = `calc. pI`,
    mw = `MW [kDa]`,
    mw_log10 = log10(mw)
  )
p2 <- ggstatsplot::ggscatterstats(
  tb_p2,
  x = pi,
  y = mw_log10,
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
ggsave(p2, file = "293T.pdf", width = 12, height = 8)
