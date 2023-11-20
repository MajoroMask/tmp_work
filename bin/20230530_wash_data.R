
library(vroom)
library(openxlsx)
library(purrr)
library(rlang)
library(stringr)
library(dplyr)
library(fs)

a <- rlang::new_environment()
a$pwd <- path("/data1/suna/work/mag0529")
setwd(a$pwd)

a$pdi <- path("./tmp/coverage/")
a$pdo <- path("./tmp/coverage_formatted") %>% dir_create()

# a$samples <- c("PPV5_E6_K", "293_PPV5_E6_K", "293_K")
a$samples <- c("293T_PPV5_E6.1")

x <- purrr::map(
  .x = a$samples,
  .f = function(smp) {
    pdi_smp <- path(a$pdi, smp)
    pdo_smp <- path(a$pdo, smp) %>% dir_create()

    # bincov
    path_bincov_in <-
      path(pdi_smp, glue::glue("{smp}_bincov.txt")) %>%
      path_real()
    path_bincov_out <-
      path(pdo_smp, glue::glue("{smp}_bincov.txt")) %>%
      file_create()
    head_bincov <-
      readr::read_lines(path_bincov_in, n_max = 3) %>%
      vroom::vroom_write_lines(file = path_bincov_out, append = FALSE)
    tb_bincov <-
      vroom::vroom(path_bincov_in, skip = 2, show_col_types = FALSE) %>%
      filter(!near(Cov, 0)) %>%
      vroom::vroom_write(file = path_bincov_out, append = TRUE)

    # covhist
    # path_covhist_in <-
    #   path(pdi_smp, glue::glue("{smp}_covhist.txt")) %>%
    #   path_real()
    # path_covhist_out <-
    #   path(pdo_smp, glue::glue("{smp}_covhist.txt")) %>%
    #   file_create()
    # tb_covhist <- vroom::vroom(path_covhist_in, show_col_types = FALSE)
    # file_copy(path_covhist_in, path_covhist_out, overwrite = TRUE)

    # covstats, rpkm & scafstats
    those <- c("covstats", "rpkm", "scafstats")
    l_tbs <-
      purrr::map(
        .x = those,
        .f = function(it) {
          pi <-
            path(pdi_smp, glue::glue("{smp}_{it}.txt")) %>%
            path_real()
          po <-
            path(pdo_smp, glue::glue("{smp}_{it}.txt")) %>%
            file_create()
          tb <-
            vroom::vroom(pi, show_col_types = FALSE, skip = if_else(it == "rpkm", 4, 0)) %>%
            filter(if_any(1, .fns = ~ !str_detect(.x, "patch of type"))) %>%
            filter(if_any(1, .fns = ~ !str_detect(.x, "unlocalized genomic scaffold"))) %>%
            filter(if_any(1, .fns = ~ !str_detect(.x, "unplaced genomic scaffold"))) %>%
            filter(if_any(1, .fns = ~ !str_detect(.x, "alternate locus group"))) %>%
            vroom::vroom_write(file = po)
          return(tb)
        }
      ) %>%
      set_names(nm = those)

    l_out <- l_tbs %>% append(values = list(covhist = tb_covhist))
    wb <- openxlsx::createWorkbook(creator = "sunamask")
    purrr::walk2(
      .x = l_out,
      .y = names(l_out),
      .f = function(x, y) {
        openxlsx::addWorksheet(wb, sheetName = y)
        openxlsx::writeData(wb, sheet = y, x = x)
        openxlsx::setColWidths(wb, sheet = y, cols = 1:ncol(x), widths = "auto")
      }
    )
    openxlsx::saveWorkbook(
      wb,
      file = path(pdo_smp, glue::glue("{smp}_formatted.xlsx")),
      overwrite = TRUE
    )
    return(l_out)
  }
)

# x[[1]]$rpkm %>% slice(1:24) %>% pull(FPKM) %>% sd()
