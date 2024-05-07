
consite <- "/data1/suna/data/tdb_playground/download"
dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

library(taxizedb)
library(rlang)
library(dplyr)

# main ----

taxizedb::tdb_cache$cache_path_set(full_path = consite)
path_sql_ncbi <- taxizedb::db_download_ncbi()

# playground ----

# link2db <- taxizedb::src_ncbi(path = path_sql_ncbi)
# taxizedb::classification(11099, db = "ncbi")
