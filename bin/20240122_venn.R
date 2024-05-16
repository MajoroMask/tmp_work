
# venn diagram

# init ----

a <- new.env(parent = emptyenv())
a$path_project <- "~/proj/tmp_work/"
a$pwd <- "/data1/suna/work/tmp_work/20240122_venn"
suppressWarnings(dir.create(a$pwd))
setwd(a$pwd)

suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(ggstatsplot))
suppressPackageStartupMessages(library(ggVennDiagram))

suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(fs))

a$pdi <- path(a$pwd, "_i/")
a$pdo <- path(a$pwd, "_o/")
dir_create(a$pdi)
dir_create(a$pdo)

# main ----

tb_input <-
  path(a$pdi, "venn_input_v6.xlsx") %>%
  openxlsx::read.xlsx(sep.names = " ") %>%
  as_tibble() %>%
  rename(acc = Accession)
l_venn <-
  map(
    .x = colnames(tb_input)[-1],
    .f = function(the_col) {
      str_output <-
        tb_input %>%
        # filter(.data[[the_col]] %in% c("High", "Peak Found")) %>%
        filter(.data[[the_col]] != 0) %>%
        pull(acc)
      return(str_output)
    }
  ) %>%
  set_names(nm = colnames(tb_input)[-1])
data_venn <- process_data(Venn(l_venn))
p_venn <-
  ggplot() +
  geom_sf(
    aes(fill = count, group = id),
    data = venn_region(data_venn),
    alpha = I(0.5)
  ) +
  geom_sf(
    aes(color = name),
    data = venn_setedge(data_venn),
    linewidth = 2,
    show.legend = TRUE
  ) +
  geom_sf_text(
    aes(label = name),
    data = venn_setlabel(data_venn)
  ) +
  geom_sf_label(
    data = venn_region(data_venn),
    mapping = aes(
      label = ifelse(
        count == 0,
        count,
        paste0(
          count, "\n",
          scales::percent(count / sum(count), accuracy = 0.1)
        )
      )
    ),
    size = 3.25
  ) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  scale_fill_gradient(low = "white", high = "#ff0000") +
  guides(
    fill = guide_colorbar(title = "prot count"),
    color = guide_legend(title = "group")
  ) +
  ggthemes::scale_color_gdocs() +
  theme_void()
p_venn

ggsave(
  p_venn, filename = path(a$pdo, "p_venn_v6.pdf"),
  width = 8, height = 6, scale = 1.5
)

# update 2024-01-23 ----

tb_input <-
  path(a$pdi, "test_input.xlsx") %>%
  openxlsx::read.xlsx(sep.names = " ") %>%
  as_tibble()
l_venn <-
  map(tb_input, na.omit) %>%
  set_names(nm = colnames(tb_input))
data_venn <- process_data(Venn(l_venn))
p_venn <-
  ggplot() +
  geom_sf(
    aes(fill = count, group = id),
    data = venn_region(data_venn),
    alpha = I(0.5)
  ) +
  geom_sf(
    aes(color = name),
    data = venn_setedge(data_venn),
    linewidth = 2,
    show.legend = TRUE
  ) +
  geom_sf_text(
    aes(label = name),
    data = venn_setlabel(data_venn)
  ) +
  geom_sf_label(
    data = venn_region(data_venn),
    mapping = aes(
      label = ifelse(
        count == 0,
        count,
        paste0(
          count, "\n",
          scales::percent(count / sum(count), accuracy = 0.1)
        )
      )
    ),
    size = 3.25
  ) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  scale_fill_gradient(low = "white", high = "#ff0000") +
  guides(
    fill = guide_colorbar(title = "prot count"),
    color = guide_legend(title = "group")
  ) +
  ggthemes::scale_color_gdocs() +
  theme_void()
p_venn
