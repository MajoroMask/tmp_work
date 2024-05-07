
# 画点图

a <- new.env(parent = emptyenv())
a$path_project <- "~/proj/tmp_work/"
a$pwd <- "/data1/suna/work/tmp_work/20240104_analysis"
suppressWarnings(dir.create(a$pwd))
setwd(a$pwd)

# suppressPackageStartupMessages(library(msa))
suppressPackageStartupMessages(library(seqinr))

suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(ggstatsplot))
suppressPackageStartupMessages(library(taxizedb))

suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(fs))

a$pdi <- path(a$pwd, "_i/")
a$pdo <- path(a$pwd, "_o/")
dir_create(a$pdi)
dir_create(a$pdo)

# utils ----

#' gsa
#'
#' group_by() %>% summarise() %>% arrange()
#'
#' @inheritParams dplyr::group_by
#' @return A tibble.
#'
#' @export
#'
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

read_xlsx <- function(...) {
  tb_out <-
    openxlsx::read.xlsx(...) %>%
    tibble::as_tibble()
  return(tb_out)
}

# functions ----

read_krakentools_combined_report <- function(input) {
  n_samples <-
    str_replace(
      readLines(input, n = 1),
      "#Number of Samples: ",
      ""
    ) %>%
    as.integer()
  tb_filename <-
    vroom::vroom(
      input,
      delim = "\t",
      col_names = c("index", "path_input"),
      col_types = "cc",
      skip = 2L,
      n_max = n_samples
    ) %>%
    mutate(
      index = str_replace(index, "^#S", ""),
      sample =
        fs::path_file(path_input) %>%
        as.character() %>%
        # str_replace("_kraken2_report.txt$", "")
        str_replace("_pe_.*$", "")
    )
  tb_data <-
    vroom::vroom(
      input,
      delim = "\t",
      col_names = c(
        "perc", "tot_all", "tot_lvl",
        paste(rep(tb_filename$sample, each = 2), c("_all", "_lvl"), sep = ""),
        "lvl_type", "taxid", "name"
      ),
      show_col_types = FALSE,
      skip = n_samples + 3L
    ) %>%
    relocate(taxid, name, .before = 1)
  return(tb_data)
}

get_taxinfo_table <- function(tax_id,
                              path_sql_ncbi,
                              path_merged_dmp,
                              reorder = FALSE) {
  requireNamespace("taxizedb", quietly = TRUE)

  # taxizedb cache setup
  assertthat::assert_that(
    fs::path_file(path_sql_ncbi) == "NCBI.sql",
    msg = "Filename of `path_sql_ncbi` must be 'NCBI.sql'!"
  )
  pd_sql_ncbi <- fs::path_dir(path_sql_ncbi)
  old_tdb_cache <- taxizedb::tdb_cache$cache_path_get()
  on.exit(taxizedb::tdb_cache$cache_path_set(full_path = old_tdb_cache))
  taxizedb::tdb_cache$cache_path_set(full_path = pd_sql_ncbi)

  # update tax_id by merged.dmp
  tb_merged_dmp <-
    path_merged_dmp %>%
    vroom::vroom(
      delim = "|",
      col_names = c("tax_id", "new", "empty"),
      col_types = "ddl",
      trim_ws = TRUE
    ) %>%
    select(-empty)
  tb_tax_id <-
    tibble(tax_id = tax_id) %>%
    left_join(tb_merged_dmp, by = "tax_id") %>%
    mutate(tax_id_updated = if_else(is.na(new), tax_id, new)) %>%
    select(tax_id, tax_id_updated)

  # get taxonomy information
  l_tax_info <-
    tb_tax_id$tax_id_updated %>%
    unique() %>%
    taxizedb::classification(db = "ncbi")
  tb_tax_info <-
    purrr::map2_dfr(
      .x = l_tax_info,
      .y = l_tax_info %>% names(),
      .f = function(tb, tax_id) {
        if (tax_id == 0) {
          tb_out <- tibble(
            tax_id = as.double(tax_id),
            name = "unclassified"
          )
        } else if (tax_id == 1) {
          tb_out <- tibble(
            tax_id = as.double(tax_id),
            name = "root"
          )
        } else if (is.null(dim(tb))) {
          tb_out <- tibble(tax_id = as.double(tax_id))
        } else {
          tb_out <- tibble(
            tax_id = tail(tb$id, 1) %>% as.double(),
            name = tail(tb$name, 1),
            rank = tail(tb$rank, 1),
            lineage_tax_id = paste(tb$id, collapse = ";"),
            lineage_name = paste(tb$name, collapse = ";"),
            lineage_rank = paste(tb$rank, collapse = ";")
          )
        }
        return(tb_out)
      }
    )
  tb_output <-
    tb_tax_id %>%
    left_join(tb_tax_info, by = c("tax_id_updated" = "tax_id"))
  return(tb_output)
}

# x         log10(all + 1)
# y         sample
# fill      level
# label     all
# facet row batch
# facet col name
plot_reads_count <- function(tb) {
  p <-
    ggplot(tb) +
    geom_bar(
      aes(y = sample, x = log10(all + 1), fill = level),
      stat = "identity"
    ) +
    geom_label(
      aes(y = sample, x = 0, label = all),
      hjust = 0, vjust = 0.5
    ) +
    ggthemes::scale_fill_gdocs(guide = "none") +
    facet_grid(
      rows = vars(batch),
      cols = vars(name),
      scales = "free",
      space = "free_y"
    ) +
    theme_bw() +
    labs(x = "# reads (log10)", y = NULL)
  return(p)
}

# x         log10(all + 1)
# y         classified_pct
# fill      level
# label     classified_pct
# facet row batch
# facet col name
plot_classified_pct <- function(tb, accu = 0.01) {
  p <-
    ggplot(tb) +
    geom_bar(
      aes(y = sample, x = classified_pct, fill = level),
      stat = "identity"
    ) +
    geom_label(
      aes(
        y = sample, x = 0,
        label = scales::percent(classified_pct, accuracy = accu)
      ),
      hjust = 0, vjust = 0.5
    ) +
    scale_x_continuous(labels = scales::percent) +
    ggthemes::scale_fill_gdocs(guide = "none") +
    facet_grid(
      rows = vars(batch),
      cols = vars(name),
      scales = "free",
      space = "free_y"
    ) +
    theme_bw() +
    labs(x = "% to classified reads", y = NULL)
  return(p)
}

# x         log10(all + 1)
# y         log10(clean_pct)
# fill      level
# label     scientific(clean_pct)
# facet row batch
# facet col name
plot_total_pct <- function(tb) {
  p <-
    ggplot(tb) +
    geom_bar(
      aes(y = sample, x = log10(clean_pct), fill = level),
      stat = "identity"
    ) +
    geom_label(
      aes(
        y = sample, x = 0.05,
        label = scales::scientific(clean_pct)
      ),
      hjust = 0, vjust = 0.5
    ) +
    scale_x_continuous(expand = expansion(add = c(0.5, 2.5))) +
    ggthemes::scale_fill_gdocs(guide = "none") +
    facet_grid(
      rows = vars(batch),
      cols = vars(name),
      scales = "free",
      space = "free_y"
    ) +
    theme_bw() +
    labs(x = "% to clean reads (log10)", y = NULL)
  return(p)
}

# x         log10(all + 1)
# y         log10(clean_pct)
# fill      level
# label     scientific(clean_pct)
# facet row batch
# facet col name
plot_umini_count <- function(tb) {
  p <-
    ggplot(tb) +
    geom_bar(
      aes(y = sample, x = log10(umini + 1), fill = level),
      stat = "identity"
    ) +
    geom_label(
      aes(y = sample, x = 0, label = umini),
      hjust = 0, vjust = 0.5
    ) +
    ggthemes::scale_fill_gdocs(guide = "none") +
    facet_grid(
      rows = vars(batch),
      cols = vars(name),
      scales = "free",
      space = "free_y"
    ) +
    theme_bw() +
    labs(x = "# unique minimizer (log10)", y = NULL)
  return(p)
}

# x         log10(all + 1)
# y         log10(clean_pct)
# fill      level
# label     scientific(clean_pct)
# facet row batch
# facet col name
plot_nmini_count <- function(tb) {
  p <-
    ggplot(tb) +
    geom_bar(
      aes(y = sample, x = log10(nmini + 1), fill = level),
      stat = "identity"
    ) +
    geom_label(
      aes(y = sample, x = 0, label = nmini),
      hjust = 0, vjust = 0.5
    ) +
    ggthemes::scale_fill_gdocs(guide = "none") +
    facet_grid(
      rows = vars(batch),
      cols = vars(name),
      scales = "free",
      space = "free_y"
    ) +
    theme_bw() +
    labs(x = "# minimizer (log10)", y = NULL)
  return(p)
}

# magic ----

a$ppl_pwd <- path("/data1/suna/work/mag_bvdv_allstar/")
a$ppl_outdir <- path(a$ppl_pwd, "outdir_test06")
a$pd_bvdv <- path(a$ppl_outdir, "Mapping/candidate_genome")

a$path_sql_ncbi <-
  path("/data1/database/taxprofiler_databases/taxdump/NCBI.sql")
taxizedb::tdb_cache$cache_path_set(full_path = fs::path_dir(a$path_sql_ncbi))

a$pd_taxdump <- path("/data1/database/taxprofiler_databases/taxdump/")

tb_samplesheet <-
  vroom(
    path(a$ppl_pwd, "samplesheet_all_20240102_add_batch.csv"),
    show_col_types = FALSE
  ) %>%
  mutate(
    sample_path = sample,
    is_ntc = str_detect(sample, "NTC"),
    level =
      case_when(
        str_detect(sample, "BVDV_\\d?E\\d_") ~
          str_replace(sample, ".*BVDV_\\d?(E\\d)_.*", "\\1"),
        str_detect(sample, "NTC") ~ "NTC",
        .default = "Non-BVDV samples"
      ) %>%
      as.factor()
  ) %>%
  arrange(batch, is_ntc, sample) %>%
  mutate(sample = fct_inorder(sample))
tb_fastp <-
  path(a$ppl_outdir, "QC_shortreads", "fastp_results.csv") %>%
  vroom::vroom(show_col_types = FALSE) %>%
  select(
    sample = sample_id,
    n_raw_reads = `summary--before_filtering.total_reads`,
    n_raw_base = `summary--before_filtering.total_bases`,
    n_clean_reads = `summary--after_filtering.total_reads`,
    n_clean_bases = `summary--after_filtering.total_bases`
  ) %>%
  mutate(
    across(
      .cols = starts_with("n_"),
      .fns = log10,
      .names = "log10_{.col}"
    )
  )

{
  tb_krakentools_combined <-
    path(a$ppl_outdir, "TaxProfile/kraken2",
         "kraken2_20231215_virus_bvdv_combined_reports.txt") %>%
    read_krakentools_combined_report() %>%
    select(
      id = taxid, name, rank = lvl_type,
      total_all = tot_all,
      total_lvl = tot_lvl
    )

  l_inputs <-
    path(a$ppl_outdir, "TaxProfile/kraken2/20231215_virus_bvdv") %>%
    dir_ls()
  l_ids <-
    l_inputs %>%
    fs::path_file() %>%
    str_replace("_pe_.*", "")
  tb_kraken2_reports <-
    map2_dfr(
      .x = l_inputs,
      .y = l_ids,
      .f = function(input_path, input_id) {
        tb_output <-
          input_path %>%
          vroom::vroom(
            delim = "\t",
            col_names = c(
              "pct", "reads_all", "reads_lvl",
              "n_mini", "n_mini_uniq", "rank", "id", "name"
            ),
            col_types = "ciiiicic",
            trim_ws = FALSE
          ) %>%
          transmute(
            id,
            name,
            rank,
            all = reads_all,
            lvl = reads_lvl,
            nmini = n_mini,
            umini = n_mini_uniq,
            sample = input_id
          )
        return(tb_output)
      }
    )
  tb_kraken2_combined <-
    tb_krakentools_combined %>%
    left_join(
      tb_kraken2_reports %>%
        pivot_wider(
          id_cols = id,
          names_from = sample,
          values_from = c(all, lvl, nmini, umini),
          names_glue = "{sample}_{.value}",
          names_vary = "slowest",
          values_fill = 0L
        ),
      by = "id"
    )
  rm(l_inputs, l_ids, tb_krakentools_combined, tb_kraken2_reports)
}

{
  tb_tmp <-
    tb_kraken2_combined %>%
    slice(1, 2) %>%
    select(name, ends_with("_all")) %>%
    select(-total_all) %>%
    pivot_longer(
      cols = -name,
      names_to = "sample",
      values_to = "reads"
    ) %>%
    pivot_wider(
      id_cols = sample,
      names_from = "name",
      values_from = "reads"
    ) %>%
    mutate(sample = str_replace(sample, "_all", "")) %>%
    rename(classified = root) %>%
    mutate(
      log10_classified = log10(classified),
      log10_unclassified = log10(unclassified)
    ) %>%
    mutate(
      pct_classified = classified / (classified + unclassified),
      pct_unclassified = unclassified / (classified + unclassified)
    )
  tb_kraken2_pct <-
    tb_samplesheet %>%
    select(sample, batch, level) %>%
    left_join(tb_tmp, by = "sample") %>%
    left_join(tb_fastp, by = "sample")
  rm(tb_tmp)

  enlist_taxids <- c(11095, 11099, 221918, 67082, 28285, 129951, 0)
  tb_tmp <-
    tb_kraken2_combined %>%
    filter(id %in% enlist_taxids) %>%
    mutate(.index = factor(id, levels = enlist_taxids)) %>%
    arrange(.index) %>%
    select(-.index) %>%
    select(name, matches("(_lvl$)|(_all$)|(_umini$)|(_nmini$)")) %>%
    select(-total_all, -total_lvl) %>%
    pivot_longer(
      cols = -name,
      names_to = "sample",
      values_to = "value"
    ) %>%
    mutate(
      type = str_replace(sample, ".*_(lvl$|all$|umini$|nmini$)", "\\1"),
      sample = str_replace(sample, "(_lvl$)|(_all$)|(_umini$)|(_nmini$)", "")
    ) %>%
    pivot_wider(
      id_cols = c(name, sample),
      names_from = "type",
      values_from = "value"
    ) %>%
    mutate(
      name =
        case_match(
          name,
          "Pestivirus" ~ "Pestivirus",
          "Bovine viral diarrhea virus 1" ~ "BVDV1",
          "Bovine viral diarrhea virus VEDEVAC" ~ "BVDV VEDEVAC",
          "BeAn 58058 virus" ~ "BAV",
          "Human adenovirus 5" ~ "HAV5",
          "Human mastadenovirus C" ~ "HMVC",
          .default = name
        ) %>%
        fct_inorder()
    )
  tb_plot_kraken2 <-
    tb_tmp %>%
    left_join(tb_kraken2_pct, by = "sample") %>%
    mutate(
      classified_pct = all / classified,
      clean_pct = all / (classified + unclassified)
    )
  rm(tb_tmp)
}

p_kraken_count_total <-
  tb_plot_kraken2 %>%
  filter(!is.na(batch)) %>%
  plot_reads_count()
p_classified_pct_total <-
  tb_plot_kraken2 %>%
  filter(!is.na(batch)) %>%
  mutate(
    classified_pct = if_else(name == "unclassified", 0, classified_pct)
  ) %>%
  plot_classified_pct()
p_clean_pct_total <-
  tb_plot_kraken2 %>%
  filter(!is.na(batch)) %>%
  plot_total_pct()
p_nmini_total <-
  tb_plot_kraken2 %>%
  filter(!is.na(batch)) %>%
  plot_nmini_count()
p_umini_total <-
  tb_plot_kraken2 %>%
  filter(!is.na(batch)) %>%
  plot_umini_count()

ggsave(
  p_kraken_count_total, filename = path(a$pdo, "p_kraken_count_0.1.pdf"),
  height = 25, width = 25
)
ggsave(
  p_classified_pct_total, filename = path(a$pdo, "p_classified_pct_0.1.pdf"),
  height = 25, width = 25
)
ggsave(
  p_clean_pct_total, filename = path(a$pdo, "p_clean_pct_0.1.pdf"),
  height = 25, width = 25
)
ggsave(
  p_nmini_total, filename = path(a$pdo, "p_nmini_total_0.1.pdf"),
  height = 25, width = 25
)
ggsave(
  p_umini_total, filename = path(a$pdo, "p_umini_total_0.1.pdf"),
  height = 25, width = 25
)

# 以上是抄的代码，现在来点奇怪的 ----

{
  enlist_taxids <- c(
    11095,
    221918,
    67082,
    1588764, 2716351, 1407493,
    2847058,
    0)
  tb_tmp <-
    tb_kraken2_combined %>%
    filter(id %in% enlist_taxids) %>%
    mutate(.index = factor(id, levels = enlist_taxids)) %>%
    arrange(.index) %>%
    select(-.index) %>%
    select(name, matches("(_lvl$)|(_all$)|(_umini$)|(_nmini$)")) %>%
    select(-total_all, -total_lvl) %>%
    pivot_longer(
      cols = -name,
      names_to = "sample",
      values_to = "value"
    ) %>%
    mutate(
      type = str_replace(sample, ".*_(lvl$|all$|umini$|nmini$)", "\\1"),
      sample = str_replace(sample, "(_lvl$)|(_all$)|(_umini$)|(_nmini$)", "")
    ) %>%
    pivot_wider(
      id_cols = c(name, sample),
      names_from = "type",
      values_from = "value"
    ) %>%
    mutate(
      name =
        case_match(
          name,
          "Pestivirus" ~ "Pestivirus",
          "Bovine viral diarrhea virus 1" ~ "BVDV1",
          "Bovine viral diarrhea virus VEDEVAC" ~ "BVDV VEDEVAC",
          "BeAn 58058 virus" ~ "BAV",
          "Human adenovirus 5" ~ "HAV5",
          "Human mastadenovirus C" ~ "HMVC",
          .default = name
        ) %>%
        fct_inorder()
    )
  tb_plot_kraken2 <-
    tb_tmp %>%
    left_join(tb_kraken2_pct, by = "sample") %>%
    mutate(
      classified_pct = all / classified,
      clean_pct = all / (classified + unclassified)
    )
  rm(tb_tmp)
}

tb_plot_kraken2 %>%
  filter(
    !is.na(batch),
    batch %in% c("240102")
  ) %>%
  plot_reads_count()
ggsave(
  path(a$pdo, "p_kraken_count_240102.pdf"),
  width = 16, height = 4
)

tb_plot_kraken2 %>%
  filter(
    !is.na(batch),
    batch %in% c("240102")
  ) %>%
  plot_nmini_count()
ggsave(
  path(a$pdo, "p_nmini_240102.pdf"),
  width = 16, height = 4
)

tb_plot_kraken2 %>%
  filter(
    !is.na(batch),
    batch %in% c("240102")
  ) %>%
  plot_umini_count()
ggsave(
  path(a$pdo, "p_umini_240102.pdf"),
  width = 16, height = 4
)
