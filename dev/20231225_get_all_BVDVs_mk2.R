
consite <- "/data1/suna/work/tmp_work/20231225_get_all_BVDVs"
dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

library(rentrez)
library(seqinr)
library(rlang)
library(glue)
library(purrr)
library(stringr)
library(tidyr)
library(dplyr)
library(fs)

library(future)
library(taxizedb)

# magic ----

a <- new_environment()

# func ----

retry <- function(expr,
                  is_error = function(x) "try-error" %in% class(x),
                  max_errors = 5L,
                  sleep = 100) {
  attempts <- 0
  try_result <- try(eval(expr))
  while (is_error(try_result)) {
    attempts <- attempts + 1
    if (attempts >= max_errors) {
      msg <- sprintf(
        "Retry: too many retries [[%s]]",
        utils::capture.output(str(try_result))
      )
      rlang::abort(msg)
    } else {
      msg <- sprintf(
        "Retry: error in attempt %i/%i [[%s]]",
        attempts,
        max_errors,
        capture.output(str(try_result))
      )
      rlang::warn(msg)
    }
    if (sleep > 0) Sys.sleep(sleep)
    try_result <- try(eval(expr))
  }
  return(try_result)
}

create_dummy_genome_info <- function(genome_info_cols = NULL) {
  genome_info_cols <- genome_info_cols %||% c(
    "uid",
    "assemblyaccession",
    "lastmajorreleaseaccession",
    "latestaccession",
    "assemblyname",
    "taxid",
    "organism",
    "assemblytype",
    "assemblystatus",
    "asmupdatedate",
    "seqreleasedate",
    "lastupdatedate",
    "releaselevel",
    "releasetype",
    "property",
    "contign50",
    "scaffoldn50",
    "ftppath_genbank",
    "ftppath_refseq"
  )

  tb_out <- purrr::map_dfc(
    .x = genome_info_cols,
    .f = ~ tibble::tibble("{.x}" := NA_character_)
  )
  return(tb_out)
}

path_ext_with_gz <- function(path) {
  first_ext <- fs::path_ext(path)
  if (first_ext == "gz") {
    second_ext <- fs::path_ext_remove(path) %>% fs::path_ext()
  } else {
    second_ext <- first_ext
    first_ext <- ""
  }
  ext_out <-
    paste(second_ext, first_ext, sep = ".") %>%
    gsub("\\.+$", "", .)
  return(ext_out)
}

download_ncbi_genome <- function(tax_id,
                                 max_errors = 5L,
                                 output_prefix,
                                 do_download = FALSE,
                                 genome_info_cols = NULL) {
  requireNamespace("rentrez", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  genome_info_cols <- genome_info_cols %||% c(
    "uid",
    "assemblyaccession",
    "lastmajorreleaseaccession",
    "latestaccession",
    "assemblyname",
    "taxid",
    "organism",
    "assemblytype",
    "assemblystatus",
    "asmupdatedate",
    "seqreleasedate",
    "lastupdatedate",
    "releaselevel",
    "releasetype",
    "property",
    "contign50",
    "scaffoldn50",
    "ftppath_genbank",
    "ftppath_refseq"
  )

  the_term <- glue::glue("txid{tax_id}[Organism:exp]")
  assemblies <- retry(
    rentrez::entrez_search(db = "assembly", term = the_term),
    max_errors = max_errors
  )

  if (assemblies$count == 0L) {
    # dummy output
    msg <- glue::glue(
      "Warning: `rentrez::entrez_search(db = 'assembly', term = {the_term})`",
      "find nothing.",
      .sep = " "
    )
    rlang::warn(msg)

    genome_filename <- NA_character_
    tb_genome_info <-
      create_dummy_genome_info(genome_info_cols) %>%
      mutate(
        tax_id = tax_id,
        ref_filename = genome_filename,
        ref_source = "NCBI_API",
        n_assemblies = 0L,
        .before = 1
      )
    return(list(genome_filename, tb_genome_info))
  }

  assemblies_summary <- retry(
    rentrez::entrez_summary(
      "assembly",
      assemblies$ids,
      always_return_list = TRUE
    ),
    max_errors = max_errors
  )
  tb_assemblies <-
    purrr::map_dfr(
      .x = assemblies_summary,
      .f = function(l) {
        l$property <-
          l$propertylist %>%
          sort() %>%
          stringr::str_c(collapse = ";")
        tb_out <- tibble::as_tibble(l[genome_info_cols])
        return(tb_out)
      }
    ) %>%
    mutate(
      is_viral_proj = case_when(
        stringr::str_detect(assemblyname, "ViralProj") ~ 1,
        TRUE ~ 99
      ),
      is_complete = case_when(
        assemblystatus == "Complete Genome" ~ 1,
        TRUE ~ 99
      ),
      across(
        .cols = c(asmupdatedate, seqreleasedate, lastupdatedate),
        .fns = lubridate::ymd_hm
      )
    ) %>%
    arrange(
      is_viral_proj,
      is_complete,
      desc(lastupdatedate),
      desc(asmupdatedate)
    )

  # download genome.fna.gz from NCBI ftp
  ftp_url <-
    rbind(tb_assemblies$ftppath_refseq, tb_assemblies$ftppath_genbank) %>%
    c() %>%
    `[`(. != "") %>%
    `[`(1) %>%
    stringr::str_c(., "/")
  ftp_contents <- retry(
    RCurl::getURL(
      url = ftp_url,
      crlf = TRUE,
      ftp.use.epsv = FALSE,
      dirlistonly = TRUE
    ),
    max_errors = max_errors
  )
  remote_filename <-
    ftp_contents %>%
    stringr::str_split(pattern = "\\n", simplify = TRUE) %>%
    as.vector() %>%
    stringr::str_subset("_genomic.fna.gz") %>%
    stringr::str_subset("_cds_from_|_rna_from_", negate = TRUE)
  tmp_xit <- path_ext_with_gz(remote_filename)
  local_filename <- glue::glue("{output_prefix}.{tmp_xit}")
  if (do_download) {
    retry(
      download.file(
        url = stringr::str_c(ftp_url, remote_filename),
        destfile = local_filename,
        quiet = TRUE
      ),
      max_errors = max_errors
    )
  }

  # prepare info table
  tb_genome_info <-
    tb_assemblies %>%
    mutate(
      .url = stringr::str_sub(ftp_url, end = -2),
      .the_row = (ftppath_refseq == .url) | (ftppath_genbank == .url)
    ) %>%
    filter(.the_row) %>%
    select(-.url, -.the_row) %>%
    mutate(
      tax_id = tax_id,
      ref_filename = local_filename,
      ref_source = "NCBI_API",
      n_assemblies = assemblies$count,
      .before = 1
    )
  return(list(local_filename, tb_genome_info))
}

taxlist2tb <- function(l) {
  purrr::map2_dfr(
    .x = l,
    .y = names(l),
    .f = ~ .x %>% mutate(parent = .y, .before = 1)
  )
}

# main ----

## get all possible BVDV list ----

x <- taxizedb::children(x = 11095) %>% taxlist2tb()
y <- taxizedb::children(x$id) %>% taxlist2tb()
z <- taxizedb::children(y$id) %>% taxlist2tb()

all_pestivirus <- c(x$id, y$id, z$id) %>% unique() %>% as.integer()

dir_create("./_o")

# if it works
# l_test <- download_ncbi_genome(
#   tax_id = 1855262,
#   output_prefix = "_o/1855262"
# )

my_future_plan <-
  future::tweak(multisession, workers = 3L) %>%
  list()
future::plan(my_future_plan)
l_results <- furrr::future_map2(
  .x = all_pestivirus,
  .y = paste0("_o/", all_pestivirus),
  .f = ~ download_ncbi_genome(
    tax_id = .x,
    output_prefix = .y,
    do_download = TRUE
  ),
  .options = furrr::furrr_options(
    packages =
      c("rentrez", "dplyr", "glue", "rlang", "stringr", "tibble",
        "lubridate", "RCurl")
  )
)
tb_combined <-
  purrr::map(
    l_results,
    .f = ~
      .x[[2]] %>%
      mutate(
        across(
          .cols = c(asmupdatedate, seqreleasedate, lastupdatedate),
          .fns = lubridate::ymd_hm
        ),
        across(
          .cols = c(contign50, scaffoldn50),
          .fns = as.integer
        )
      )
  ) %>%
  bind_rows()
tb_output_extended <-
  tb_combined %>%
  filter(!is.na(ref_filename)) %>%
  group_by(uid) %>%
  group_map(
    .f = ~ slice_head(.x, n = 1L),
    .keep = TRUE
  ) %>%
  bind_rows() %>%
  mutate(
    name = str_replace(organism, fixed(" (viruses)"), ""),
    path = fs::path_abs(ref_filename)
  )

# 重做！

l_results <- furrr::future_map2(
  .x = tb_output_extended$tax_id,
  .y = glue::glue("_o/genome_{tb_output_extended$tax_id}"),
  .f = ~ download_ncbi_genome(
    tax_id = .x,
    output_prefix = .y,
    do_download = TRUE
  ),
  .options = furrr::furrr_options(
    packages =
      c("rentrez", "dplyr", "glue", "rlang", "stringr", "tibble",
        "lubridate", "RCurl")
  )
)

tb_combined <-
  purrr::map(
    l_results,
    .f = ~
      .x[[2]] %>%
      mutate(
        across(
          .cols = c(asmupdatedate, seqreleasedate, lastupdatedate),
          .fns = lubridate::ymd_hm
        ),
        across(
          .cols = c(contign50, scaffoldn50),
          .fns = as.integer
        )
      )
  ) %>%
  bind_rows()
tb_output_extended <-
  tb_combined %>%
  filter(!is.na(ref_filename)) %>%
  group_by(uid) %>%
  group_map(
    .f = ~ slice_head(.x, n = 1L),
    .keep = TRUE
  ) %>%
  bind_rows() %>%
  mutate(
    name = str_replace(organism, fixed(" (viruses)"), ""),
    path = fs::path_abs(ref_filename)
  )

vroom::vroom_write(
  tb_output_extended,
  file = path("_o/refsheet_BVDV_allstar_extendend.csv"),
  delim = ","
)
vroom::vroom_write(
  tb_output_extended %>% select(tax_id, reference = path),
  file = path("_o/refsheet_BVDV_allstar.csv"),
  delim = ","
)
