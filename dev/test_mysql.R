
library(DBI)
library(RMySQL)
library(dbplyr)
library(dplyr)

con <- DBI::dbConnect(
  RMySQL::MySQL(),
  host = "10.10.110.201",
  port = 3306,
  dbname = "thatdb",
  user = "suna",
  password = "4561"
)

tb_ngs_sample <-
  tbl(con, "ams_ngs_sample") %>%
  collect()

dbDisconnect(con)
