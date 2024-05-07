
consite <- "/data1/suna/work/tmp_work/20231205_get_all_BVDVs"
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

# magic ----

a <- new_environment()
a$enlist_taxids <- c(
  11099, 11100, 31656, 125217, 140298, 140303, 144353, 144355, 217536, 229907,
  268305, 981676, 1911095,

  64004, 64005, 64006, 64007, 68642, 68643, 68644, 68645, 68646, 68647, 68648,
  68649, 68650, 68653, 68654, 68655, 68656, 68657, 68658, 68659, 69043, 69044,
  82470, 83397, 99470, 99471, 99472, 99474, 99475, 99476, 99477, 99542, 99543,
  116076, 116077, 116078, 116079, 116080, 116081, 116082, 116083, 120950,
  121430, 144356, 144357, 144358, 144359, 144360, 144361, 144362, 144363,
  144364, 144365, 144366, 144367, 144368, 144369, 144370, 144371, 145170,
  145171, 145172, 145173, 145174, 145175, 145176, 145177, 145178, 145179,
  145180, 145181, 145182, 145183, 145184, 145185, 145186, 145187, 145188,
  145189, 145190, 145191, 145192, 145193, 145194, 145195, 145196, 145197,
  145198, 145199, 145200, 145201, 145202, 145203, 145204, 145205, 145206,
  145207, 145208, 145209, 145210, 145211, 145212, 145213, 145214, 145215,
  145216, 145217, 145218, 145219, 145220, 145221, 145222, 145223, 145224,
  145225, 145226, 145227, 145228, 145229, 158474, 158476, 158477, 158901,
  158902, 158903, 158904, 158905, 158906, 158907, 158908, 158909, 158910,
  158911, 158912, 158913, 158914, 158915, 158916, 158917, 158918, 158919,
  158920, 158921, 158922, 158923, 158924, 158925, 217899, 217900, 217901,
  217902, 217903, 217904, 217905, 217906, 217907, 217908, 217909, 217910,
  217911, 217912, 217913, 217914, 217915, 217916, 217917, 217918, 217919,
  217920, 217921, 217922, 217923, 217924, 217925, 217926, 217927, 217928,
  217929, 217930, 217931, 217932, 217933, 217934, 217935, 221916, 221917,
  221918, 230462, 230463, 230464, 233883, 233884, 233885, 233886, 233887,
  233888, 233889, 258254, 262958, 262959, 262960, 262961, 262962, 262963,
  312033, 312034, 312035, 365500, 365501, 365502, 365503, 365504, 365505,
  365506, 365507, 365508, 365509, 365510, 365511, 365512, 365513, 365514,
  664649, 664654, 1855262, 1855263, 1855264, 1855265, 1855266, 1855267,
  1855268, 1855269, 2590438
)

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
  retry(
    download.file(
      url = stringr::str_c(ftp_url, remote_filename),
      destfile = local_filename,
      quiet = TRUE
    ),
    max_errors = max_errors
  )

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

# main ----

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
  .x = a$enlist_taxids,
  .y = paste0("_o/", a$enlist_taxids),
  .f = ~ download_ncbi_genome(
    tax_id = .x,
    output_prefix = .y
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
tb_output <-
  tb_combined %>%
  filter(!is.na(ref_filename)) %>%
  group_by(uid) %>%
  group_map(
    .f = ~ slice_head(.x, n = 1L),
    .keep = TRUE
  ) %>%
  bind_rows() %>%
  transmute(
    name = str_replace(organism, fixed(" (viruses)"), ""),
    path = fs::path_abs(ref_filename)
  )
vroom::vroom_write(
  tb_output,
  file = path("_o/refsheet_BVDV_allstar.csv"),
  delim = ","
)
