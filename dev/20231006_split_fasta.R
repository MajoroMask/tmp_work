#!/usr/local/bin/R

# Read fasta from kraken2 library.fna
# split by genome
# write to separate .fa.gz files

library(rlang)
library(stringr)
library(purrr)
library(seqinr)
library(fs)

wd <- "/data1/suna/data/tdb_playground/kmcp"
setwd(wd)

seq_all <- read.fasta("library.fna", seqtype = "DNA", as.string = TRUE)
tax_id <-
  getName(seq_all) %>%
  stringr::str_replace(".*\\|(\\d+)\\|.*", "\\1")
tax_id_uniq <- unique(tax_id)
gz_filenames <-
l_seq_out <- walk(
  .x = tax_id_uniq,
  .f = function(the_id) {
    l_sub <- seq_all[[tax_id == the_id]]
    seqinr::write.fasta(
      # TODO
      # see https://bioinf.shenwei.me/kmcp/database/#building-custom-databases
    )
  },
  seq_all = seq_all,
  tax_id = tax_id
)
