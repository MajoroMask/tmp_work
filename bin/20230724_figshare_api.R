
# Download data from figshare API
# All available data from this artical:
# https://www.nature.com/articles/s41592-021-01141-3

consite <- "/data1/suna/proj/nmethod_benchmarking_2021"
R_proj <- "~/proj/tmp_work"
# dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

renv::activate(R_proj)

library(httr)
library(jsonlite)

library(future)
library(furrr)

library(glue)
library(rlang)
library(purrr)
library(dplyr)
library(fs)

# magic ----

proj_id <- "79916"
pdo <- "./fastq" %>% dir_create()

# main ----

proj_articles <-
  glue("https://api.figshare.com/v2/projects/{proj_id}/articles") %>%
  httr::GET() %>%
  httr::content()
articles_api <-
  proj_articles %>%
  purrr::map_chr(.f = ~ .$url_public_api)
tb_files <-
  articles_api %>%
  purrr::map_dfr(
    .f = ~
      glue("{.x}/files") %>%
      httr::GET() %>%
      httr::content() %>%
      purrr::map_dfr(.f = ~ as_tibble(.))
  ) %>%
  mutate(
    destfile = path(pdo, name) %>% path_rel()
  )

future::plan(multisession, workers = 8)
furrr::future_walk2(
  .x = tb_files$id,
  .y = tb_files$destfile,
  .f = ~
    httr::GET(
      url = glue("https://api.figshare.com/v2/file/download/{.x}"),
      write_disk(.y, overwrite = TRUE)
    )
)
