
# 把20231206_analysis.Rmd里的代码整理一下，画图贴到PPT里

# init ----

a <- new.env(parent = emptyenv())
a$path_project <- "~/proj/tmp_work/"
a$pwd <- "/data1/suna/work/tmp_work/20231227_analysis"
suppressWarnings(dir.create(a$pwd))
setwd(a$pwd)

suppressPackageStartupMessages(library(msa))
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

# magic ----

a$ppl_pwd <- path("/data1/suna/work/mag1205/")
a$ppl_outdir <- path(a$ppl_pwd, "outdir_test03")
a$pd_bvdv <- path("/data1/suna/work/tmp_work/20231205_get_all_BVDVs/_o")

a$path_sql_ncbi <-
  path("/data1/database/taxprofiler_databases/taxdump/NCBI.sql")
taxizedb::tdb_cache$cache_path_set(full_path = fs::path_dir(a$path_sql_ncbi))

a$pd_taxdump <- path("/data1/database/taxprofiler_databases/taxdump/")

# fonts

a$font_family <- "sarasa-term-sc-nerd-regular"
a$font_regular <-
  path(
    "/etc", "rstudio", "fonts", "sarasa\ term\ sc\ nerd", "400",
    "sarasa-term-sc-nerd-regular.ttf"
  ) %>%
  path_real()
sysfonts::font_add(
  family = a$font_family,
  regular = a$font_regular
)
showtext::showtext_auto()

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

# function ----

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

# 数据 ----

## 0825/1202 ----

# samplesheet
tb_samplesheet <-
  vroom(
    path(a$ppl_pwd, "samplesheet_extended.csv"),
    show_col_types = FALSE
  ) %>%
  mutate(
    sample_path = sample,
    sample = str_replace(sample_path, "H7N9_E9_BVDV_7", "BVDV_"),
    batch = str_replace(sample, ".*_", ""),
    is_ntc = str_detect(sample, "NTC"),
    level =
      str_replace(sample, "BVDV_(.*?)_.*", "\\1") %>%
      str_replace("NTC.*", "NTC") %>%
      fct_inorder()
  ) %>%
  arrange(batch, is_ntc, sample) %>%
  mutate(sample = fct_inorder(sample))

# fastp
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
  mutate(sample = str_replace(sample, "H7N9_E9_BVDV_7", "BVDV_"))

# kraken2

path_krakentool_combined_report <-
  path(a$ppl_outdir, "Taxonomy", "kraken2_combined_reports.txt")
tb_kraken_input <-
  path_krakentool_combined_report %>%
  read_krakentools_combined_report()
tb_tax_anno_kraken <- get_taxinfo_table(
  tax_id = tb_kraken_input$taxid,
  path_sql_ncbi = a$path_sql_ncbi,
  path_merged_dmp = fs::path(a$pd_taxdump, "merged.dmp"),
  reorder = FALSE
)
tb_kraken_final <-
  tb_kraken_input %>%
  rename(name_ori = name) %>%
  left_join(tb_tax_anno_kraken, by = c("taxid" = "tax_id")) %>%
  relocate(
    tax_id_updated,
    .after = taxid
  ) %>%
  relocate(
    name, rank, lineage_tax_id, lineage_name, lineage_rank,
    .after = name_ori
  )
tb_tmp <-
  tb_kraken_final %>%
  slice(1, 2) %>%
  select(name, ends_with("_all")) %>%
  select(-tot_all) %>%
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
  mutate(
    sample =
      str_replace(sample, "_all", "") %>%
      str_replace("H7N9_E9_BVDV_7", "BVDV_"),
  ) %>%
  rename(classified = root) %>%
  mutate(pct = classified / (classified + unclassified)) %>%
  mutate(
    log10_classified = log10(classified),
    log10_unclassified = log10(unclassified)
  )
tb_kraken_pct <-
  tb_samplesheet %>%
  select(sample, batch, level) %>%
  left_join(tb_tmp, by = "sample") %>%
  left_join(tb_fastp, by = "sample")

# kraken2 subset
tb_tmp <-
  tb_kraken_final %>%
  filter(taxid %in% c(11095, 11099, 67082, 28285)) %>%
  select(name, ends_with("_all")) %>%
  select(-tot_all) %>%
  pivot_longer(
    cols = -name,
    names_to = "sample",
    values_to = "reads"
  ) %>%
  mutate(
    sample =
      str_replace(sample, "_all", "") %>%
      str_replace("H7N9_E9_BVDV_7", "BVDV_"),
  ) %>%
  mutate(
    name =
      case_match(
        name,
        "Bovine viral diarrhea virus 1" ~ "BVDV1",
        "BeAn 58058 virus" ~ "BAV",
        "Human adenovirus 5" ~ "HAV5",
        .default = name
      ) %>%
      fct_inorder()
  )
tb_plot_pestivirus <-
  tb_tmp %>%
  left_join(tb_kraken_pct, by = "sample") %>%
  mutate(
    classified_pct = reads / classified,
    clean_pct = reads / n_clean_reads
  )
