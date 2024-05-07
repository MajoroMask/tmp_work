
# init ----

library(rlang)
library(vroom)
library(purrr)
library(tidyr)
library(stringr)
library(forcats)
library(dplyr)
library(ggplot2)
library(fs)

setwd("~/work/tmp_work/20230506_batch_comare/")

a <- rlang::new_environment()
a$pdi <- path("./_i/")
a$pdo <- path("./_o/")

# function ----

read_fastqc_table <- function(tb) {
  batch <- tb$batch %>% unique() %>% sort()
  tb_fastqc <- purrr::map_dfr(
    .x = batch,
    .f = ~
      fs::path(a$pdi, .x, "fastqc_results.csv") %>%
      vroom(show_col_types = FALSE)
  )
  tb_out <-
    tb %>%
    left_join(
      y = tb_fastqc,
      by = c("sample" = "sample_id")
    )
  return(tb_out)
}

read_bracken_output <- function(tb) {
  l <-
    pmap(
      .l = tb %>% select(batch, sample),
      .f = function(batch, sample) {

        if (sample == "lib102") {
          browser()
        }

        tb_out <-
          glue::glue("./_i/{batch}/bracken/{sample}/{sample}_S.bracken_output.txt") %>%
          fs::path() %>%
          vroom::vroom(show_col_types = FALSE) %>%
          select(
            name, taxonomy_id, taxonomy_lvl,
            new_est_reads, fraction_total_reads
          ) %>%
          rename(
            "{sample}_num" := new_est_reads,
            "{sample}_frac" := fraction_total_reads,
          )
        return(tb_out)
      }
    )
  tb_out <-
    purrr::reduce(
      .x = l,
      .f = full_join,
      by = c("name", "taxonomy_id", "taxonomy_lvl")
    ) %>%
    mutate(
      across(
        .col = -c(name, taxonomy_id, taxonomy_lvl),
        .fns = ~ tidyr::replace_na(.x, replace = 0)
      )
    )
  return(tb_out)
}

# input ----

tb_samples <- path(a$pdi, "samples.tsv") %>% vroom(show_col_types = FALSE)
tb_qc <- read_fastqc_table(tb_samples)

## base/reads count

# reads数也差太多了吧……

tb_p1 <-
  tb_qc %>%
  select(
    sample, seq_id, host, batch,
    matches("total_sequence.*trimmed", perl = TRUE)
  ) %>%
  rename_with(
    .cols = matches("total_sequence.*trimmed", perl = TRUE),
    .fn = ~ str_replace(.x, ".*(read[12]).*", "\\1")
  ) %>%
  pivot_longer(
    cols = c(read1, read2),
    names_to = "read",
    values_to = "value"
  ) %>%
  mutate(
    sample = fct_inorder(sample),
    host = fct_inorder(host)
  )
p1_log10 <-
  ggplot(tb_p1) +
  geom_bar(
    aes(sample, value, fill = host),
    stat = "identity"
  ) +
  scale_y_log10(n.breaks = 10) +
  ggthemes::scale_fill_calc() +
  facet_grid(
    rows = vars(read),
    cols = vars(batch),
    scales = "free_x",
    space = "free"
  ) +
  ggthemes::theme_calc() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
p1_regular <-
  ggplot(tb_p1) +
  geom_bar(
    aes(sample, value, fill = host),
    stat = "identity"
  ) +
  scale_y_continuous(n.breaks = 10) +
  ggthemes::scale_fill_calc() +
  facet_grid(
    rows = vars(read),
    cols = vars(batch),
    scales = "free_x",
    space = "free"
  ) +
  ggthemes::theme_calc() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# bracken结果整理 ----

tb_bracken_ori <- read_bracken_output(tb_samples)

tb_bracken_long <-
  tb_bracken_ori %>%
  pivot_longer(
    cols = -c(name, taxonomy_id, taxonomy_lvl),
    names_pattern = "(.+)_(num|frac)",
    names_to = c("sample", "data_type")
  )
x <-
  tb_bracken_long %>%
  filter(
    str_detect(name, regex("ovine respirovirus", ignore_case = TRUE))
  )
x %>% filter(data_type == "num") %>% arrange(desc(value)) %>% left_join(tb_samples, by = "sample") %>% View()
y <-
  tb_bracken_long %>%
  filter(
    str_detect(name, regex("ovine m", ignore_case = TRUE))
  )
y %>% filter(data_type == "num") %>% arrange(desc(value)) %>% left_join(tb_samples, by = "sample") %>% View()
# centrifuge ----

read_centrifuge_output <- function(tb) {
  l <-
    pmap(
      .l = tb %>% select(batch, sample),
      .f = function(batch, sample) {

        if (sample == "lib102") browser()

        tb_out <-
          glue::glue("./_i/{batch}/centrifuge/{sample}/report.txt") %>%
          fs::path() %>%
          vroom::vroom(show_col_types = FALSE) %>%
          filter(taxRank == "species") %>%
          select(
            name,
            taxonomy_id = taxID,
            taxonomy_lvl = taxRank,
            n_reads = numReads,
            n_reads_uniq = numUniqueReads,
            abundance
          ) %>%
          rename(
            "{sample}_num" := n_reads,
            "{sample}_num_uniq" := n_reads_uniq,
            "{sample}_frac" := abundance
          )
        return(tb_out)
      }
    )
}

l_centrifuge <- read_centrifuge_output(tb_samples)

tb_centrifuge_homo <-
  map_dfr(
    l_centrifuge,
    .f = ~
      .x %>%
      filter(taxonomy_id  == 9606) %>%
      rename_with(
        .cols = matches("_(num|frac)", perl = TRUE),
        .fn = ~ str_replace(.x, ".*_(num|frac)", "\\1")
      ) %>%
      select(name, num, frac)
  ) %>%
  bind_cols(tb_samples)
tb_p_centri_homo <-
  tb_centrifuge_homo %>%
  pivot_longer(
    cols = c(num, frac),
    names_to = "data_type",
    values_to = "value"
  ) %>%
  mutate(
    sample = fct_inorder(sample),
    host = fct_inorder(host)
  )
p_centrifuge_homo_num <-
  tb_p_centri_homo %>%
  filter(data_type == "num") %>%
  ggplot() +
  geom_bar(
    aes(sample, value, fill = host),
    stat = "identity",
    alpha = I(0.9)
  ) +
  scale_y_log10(n.breaks = 10) +
  ggthemes::scale_fill_calc() +
  facet_grid(
    cols = vars(batch),
    scales = "free",
    space = "free"
  ) +
  ggthemes::theme_calc() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(y = "n reads")
p_centrifuge_homo_num
p_centrifuge_homo_frac <-
  tb_p_centri_homo %>%
  filter(data_type == "frac") %>%
  ggplot() +
  geom_bar(
    aes(sample, value, fill = host),
    stat = "identity",
    alpha = I(0.9)
  ) +
  scale_y_continuous(labels = scales::percent) +
  ggthemes::scale_fill_calc() +
  facet_grid(
    cols = vars(batch),
    scales = "free",
    space = "free"
  ) +
  ggthemes::theme_calc() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(y = "reads fraction")
p_centrifuge_homo_frac

#
