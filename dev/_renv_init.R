
# renv::deactivate()
chooseCRANmirror()  # use bfsu mirror
chooseBioCmirror()  # use TUNA mirror

# .libPaths()
# renv:::renv_paths_root()
# BiocManager::repositories()

renv::init(
  repos = BiocManager::repositories(),
  bioconductor = TRUE
)

renv::update(packages = "renv")

renv::install("RSQLite")
renv::install("tidyverse")
renv::install("tidymodels")
