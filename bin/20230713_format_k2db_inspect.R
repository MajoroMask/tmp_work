
consite <- "/data1/suna/work/tmp_work/20230713_k2db_inspect"
# dir.create(consite, recursive = TRUE, mode = "0775")
setwd(consite)

library(rlang)
library(vroom)
library(purrr)
library(stringr)
library(tidyr)
library(dplyr)
library(fs)

# utils ----

gsa <- function(.data, ..., .sort = TRUE) {
  requireNamespace("dplyr", quietly = TRUE)
  .data %>%
    group_by(...) %>%
    tally(sort = .sort) %>%
    ungroup() %>%
    # mutate(pct = round(n / sum(n), digits = 3))
    mutate(pct = scales::percent(round(n / sum(n), digits = 3)))
}

# function ----

xit_fun <- function(tb, current_level, id, accumulate_var, accumulate_full) {
  upper_level <- current_level - 2
  the_name <- glue::glue("{accumulate_var}_{current_level}")

  tb_sub <-
    tb %>%
    filter(cascade_level %in% c(current_level, upper_level)) %>%
    transmute(
      id = .data[[id]],
      accu_var = .data[[accumulate_var]],
      accu_full = .data[[accumulate_full]],
      cascade_level = cascade_level,
      "{the_name}" := if_else(
        cascade_level == current_level,
        accu_var,
        NA_character_
      ),
      accu_full_new = if_else(
        cascade_level == current_level,
        str_c(accu_full, accu_var, sep = ";"),
        accu_full
      )
    ) %>%
    tidyr::fill({{the_name}}) %>%
    filter(cascade_level == current_level) %>%
    select(id, all_of(the_name), accu_full_new)

  tb_new <-
    tb %>%
    transmute(
      id = .data[[id]],
      accu_var = .data[[accumulate_var]],
      accu_full = .data[[accumulate_full]],
      cascade_level = cascade_level
    ) %>%
    left_join(tb_sub, by = "id") %>%
    mutate(
      accu_full = case_when(
        cascade_level == current_level ~ accu_full_new,
        cascade_level < current_level ~ accu_full,
        cascade_level > current_level ~ NA_character_
      )
    ) %>%
    select(-accu_full_new) %>%
    tidyr::fill({{the_name}}, accu_full) %>%
    mutate(
      across(
        .cols = all_of(the_name),
        .fns = ~ if_else(cascade_level < current_level, NA_character_, .x)
      )
    )
  tb_out <-
    tb %>%
    mutate(
      "{accumulate_full}" := tb_new %>% pull(accu_full),
      "{the_name}" := tb_new %>% pull({{the_name}})
    )
  return(tb_out)
}

# magic ----

k2_inspect_colnames <-
  c(
    "frac_clade", "mini_clade", "mini_direct",
    "tax_level", "tax_id", "tax_info"
  )
fct_level_tax <-
  c(
    "D", "D1", "D2", "D3", "D4", "D5",
    "K",
    "P", "P1",
    "C", "C1", "C2",
    "O", "O1",
    "F", "F1", "F2",
    "G", "G1", "G2",
    "S", "S1", "S2", "S3"
  )

# main ----

tb_input <-
  path("inspect.txt") %>%
  vroom(
    delim = "\t",
    skip = 7L,
    col_names = k2_inspect_colnames,
    trim_ws = FALSE,
    show_col_types = FALSE
  )
tb_format <-
  tb_input %>%
  select(tax_info, tax_level) %>%
  mutate(
    tax_name = stringr::str_trim(tax_info),
    tax_name_full = "root",
    tax_name_0 = "root",
    tax_level_full = "R",
    tax_level_0 = "R",
    cascade_level =
      stringr::str_replace_all(tax_info, "^(\\s*).*", "\\1") %>%
      stringr::str_count()
  ) %>%
  filter(tax_name != "environmental samples") %>%
  select(
    cascade_level,
    tax_name, tax_name_full, tax_name_0,
    tax_level, tax_level_full, tax_level_0
  )
cascade_levels <-
  tb_format %>%
  pull(cascade_level) %>%
  unique() %>%
  sort() %>%
  tail(-1)

tb_tax_name <- purrr::reduce(
  .x = cascade_levels,
  .f = xit_fun,
  .init = tb_format,
  id = "tax_name",
  accumulate_var = "tax_name",
  accumulate_full = "tax_name_full"
)
tb_tax_level <- purrr::reduce(
  .x = cascade_levels,
  .f = xit_fun,
  .init = tb_format,
  id = "tax_name",
  accumulate_var = "tax_level",
  accumulate_full = "tax_level_full"
)
tb_output <-
  tb_input %>%
  mutate(tax_name = stringr::str_trim(tax_info)) %>%
  left_join(
    tb_tax_name %>% select(tax_name, starts_with("tax_name_")),
    by = "tax_name"
  ) %>%
  left_join(
    tb_tax_level %>% select(tax_name, starts_with("tax_level_")),
    by = "tax_name"
  ) %>%
  select(-tax_name) %>%
  relocate(tax_level_full, .after = tax_name_full)
vroom::vroom_write(
  tb_output,
  file = "inspect_flattened.tsv",
  na = "-"
)
