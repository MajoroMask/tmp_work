
consite <- "/data1/suna/work/tmp_work/20230912_download_s3"
dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

R_proj <- "~/proj/tmp_work"
renv::activate(R_proj)

library(aws.s3)
library(rlang)
library(dplyr)
library(fs)

a <- new_environment()
a$pdi <- path("_i") %>% dir_create()
a$pdo <- path("_o") %>% dir_create()

# get this from https://registry.opendata.aws/aws-igenomes/
a$bucket <- "s3://ngi-igenomes/"
a$region <- "eu-west-1"

bucket_exists(a$bucket, region = a$region)
tb_bucket_files <-
  get_bucket_df(a$bucket, region = a$region) %>%
  as_tibble()
save_object(
  "README.md", bucket = a$bucket, region = a$region,
  file = path(a$pdo, "README.md")
)
save_object(
  "test-data/taxprofiler/db_mOTU.tar.gz",
  bucket = a$bucket, region = a$region,
  file = path(a$pdo, "db_mOTU.tar.gz")
)

# TODO 下载mOTU的test reference
