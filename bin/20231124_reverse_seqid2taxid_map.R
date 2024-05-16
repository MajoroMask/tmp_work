
consite <- "/data1/suna/work/tmp_work/20231124_reverse_seqid2taxid_map"
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

catch_error <- function(expr,
                        is_error = function(x) "try-error" %in% class(x)) {
  try_result <- try(eval(expr))
  while (is_error(try_result)) {
    browser()
  }
  try_result
}

pipe_end <- function(x) x

gsa <- function(.data, ..., .sort = TRUE) {
  requireNamespace("dplyr", quietly = TRUE)
  .data %>%
    group_by(...) %>%
    tally(sort = .sort) %>%
    ungroup() %>%
    # mutate(pct = round(n / sum(n), digits = 3))
    mutate(pct = scales::percent(round(n / sum(n), digits = 3)))
}

# main ----

arg <- rlang::new_environment()
arg$path_id_map <-
  path("/data1/database/taxprofiler_databases/kraken2/20230926_virus/seqid2taxid.map")

tb_idmap <-
  vroom::vroom(
    arg$path_id_map,
    delim = "\t",
    col_names = c("seq_id", "tax_id"),
    col_types = "ci"
  ) %>%
  mutate(
    nucl_id = str_replace(seq_id, ".*\\|", ""),
    .after = seq_id
  )

# Part1：tax_id和seq_id一一对应
tax_with_unique_assembly <-
  tb_idmap %>%
  group_by(tax_id) %>%
  summarise(
    group_size = n()
  ) %>%
  ungroup() %>%
  filter(group_size == 1L) %>%
  pull(tax_id)
tb_tax2seq_uniq <-
  tb_idmap %>%
  filter(tax_id %in% tax_with_unique_assembly) %>%
  select(tax_id, seq_id, nucl_id)

# Part 2：一个tax_id下有多于一个seq_id的情况
#   通过assembly归拢seq_id，先通过seq_id获取assembly_id，
#   再通过assembly_id归拢隶属同一个assembly的多个seq_id
tb_idmap_dup <-
  tb_idmap %>%
  filter(!tax_id %in% tax_with_unique_assembly)
# 使用NCBI API可能会抽疯，所以用这种分批次的方式
l_links_all <- purrr::map(
  .x = tb_idmap_dup$nucl_id[1:100],
  .f = link_nucleotide2get_assembly
)
for (i in seq(101, nrow(tb_idmap_dup), 50)) {
  l_links_all_tmp <- purrr::map(
    .x = tb_idmap_dup$nucl_id[i:pmin(i + 49, nrow(tb_idmap_dup))],
    .f = link_nucleotide2get_assembly
  )
  print(i)
  l_links_all <- c(l_links_all, l_links_all_tmp)
}

# 发现了一个nucleotide属于多个assemblies的情况，
# 看了一下，是同一个assembly有多个版本
tb_asmbl_per_nucl <- map2_dfr(
  .x = tb_idmap_dup$nucl_id,
  .y = l_links_all,
  .f = function(nucl_id, elink) {
    tibble(
      nucl_id = nucl_id,
      n_asmbl = length(elink$links$nuccore_assembly),
      asmbl_id_str = paste(elink$links$nuccore_assembly, collapse = ";")
    )
  }
)

# Part 3：通过nucleotide搭桥，筛选出一个tax_id对应多个assembly_id的情况。
tb_idmap_dup_ext <-
  tb_idmap_dup %>%
  left_join(tb_asmbl_per_nucl, by = "nucl_id")
tb_idmap_dup_ext_dedup <-
  tb_idmap_dup_ext %>%
  group_by(tax_id) %>%
  summarise(
    xit = length(unique(asmbl_id_str)) == 1,
    asmbl_id_str_regrp =
      asmbl_id_str %>%
      str_split(pattern = ";") %>%
      unlist() %>%
      unique() %>%
      sort() %>%
      paste(collapse = ";")
  ) %>%
  ungroup() %>%
  filter(!xit) %>%
  rowwise() %>%
  group_map(
    .f = ~ tibble(
      tax_id = .x$tax_id,
      asmbl_id = str_split(.x$asmbl_id_str_regrp, pattern = ";") %>% unlist()
    )
  ) %>%
  bind_rows()

# 获取这些assemblies的summary信息，优先级排序
l_asmbl_summary <-
  purrr::map(
    .x = tb_idmap_dup_ext_dedup$asmbl_id,
    .f = ~ retry(entrez_summary(db = "assembly", id = .x))
  )
enlist_cols <-
  purrr::map(
    .x = l_asmbl_summary,
    .f = ~ map_int(.x, .f = length) %>% `[`(. == 1L) %>% names()
  ) %>%
  purrr::reduce(.f = intersect)
tb_asmbl_summary <-
  purrr::map_dfr(
    l_asmbl_summary,
    .f = ~ .x[enlist_cols] %>% as_tibble()
  ) %>%
  left_join(tb_idmap_dup_ext_dedup, by = c("uid" = "asmbl_id")) %>%
  relocate(tax_id, .before = 1) %>%
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

# 为一个tax_id指定一个优先选择的assembly_id
tb_asmbl_summary_enlist <-
  tb_asmbl_summary %>%
  group_by(tax_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(tax_id, asmbl_id = uid) %>%
  left_join(
    tb_asmbl_per_nucl %>%
      rowwise() %>%
      group_map(
        .f = ~ tibble(
          asmbl_id = str_split(.$asmbl_id_str, pattern = ";") %>% unlist(),
          nucl_id = .$nucl_id
        )
      ) %>%
      bind_rows() %>%
      distinct(),
    by = "asmbl_id"
  )

# 最终的结果由三部分组成：
# - 在idmap中，tax_id和seq_id一一对应，可以直接反向索引的部分
# - tax_id在tb_idmap_dup_ext$tax_id中，但不在tb_idmap_dup_ext_dedup$tax_id中
#   这部分是一个tax_id只有一个assembly，但包含多个seq_id的情况，
#   即不是 complete genome
# - 在tb_asmbl_summary_enlist中的，这部分是一个tax_id对应多个assemblies，
#   通过优先级选择一个assembly。
tb_final_xits <-
  bind_rows(
    tb_tax2seq_uniq %>% select(tax_id, seq_id),
    tb_idmap %>%
      filter(
        tax_id %in% tb_idmap_dup_ext$tax_id,
        !(tax_id %in% tb_idmap_dup_ext_dedup$tax_id)
      ) %>%
      select(tax_id, seq_id),
    tb_asmbl_summary_enlist %>%
      left_join(
        tb_idmap %>% select(nucl_id, seq_id),
        by = "nucl_id"
      ) %>%
      select(tax_id, seq_id)
  )
vroom::vroom_write(
  tb_final_xits,
  file = "taxid2seqid.map",
  delim = "\t",
  col_names = FALSE
)
