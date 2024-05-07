
# venn diagram, party of two

# init ----

a <- new.env(parent = emptyenv())
a$path_project <- "~/proj/tmp_work/"
a$pwd <- "/data1/suna/work/tmp_work/20240125_analysis"
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
  path(a$pdi, "20240104-JM(1).xlsx") %>%
  openxlsx::read.xlsx(sep.names = " ", na.strings = c("NA", "")) %>%
  as_tibble() %>%
  rename(acc = Accession)

# venn ----

l_venn <-
  map(
    .x = colnames(tb_input)[-1],
    .f = function(the_col) {
      str_output <-
        tb_input %>%
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

# ggsave(
#   p_venn, filename = path(a$pdo, "p_venn_v6.pdf"),
#   width = 8, height = 6, scale = 1.5
# )

# cor ----

m_cor_input <-
  tb_input %>%
  select(-`Ratio:JS005-HCP/JS005-FS`) %>%
  na.omit() %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "acc") %>%
  as.matrix()
cor(m_cor_input)

tb_p <-
  tb_input %>%
  rename(
    ratio = `Ratio:JS005-HCP/JS005-FS`,
    abun1 = `JS005-FS`,
    abun2 = `JS005-HCP`
  ) %>%
  filter(
    !is.na(abun1),
    !is.na(abun2)
  ) %>%
  mutate(
    color =
      case_when(
        log2(ratio) > 2 ~ "up",
        log2(ratio) < -2 ~ "down",
        TRUE ~ "not sign"
      ) %>%
      factor(levels = c("not sign", "down", "up"))
  ) %>%
  arrange(color) %>%
  mutate(color = factor(color, levels = c("up", "not sign", "down")))

major_breaks_at <- c(
  0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000,
  1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8
)
minor_breaks_at <- c(
  25, 75, 250, 750, 2500, 7500, 25000, 75000,
  2.5e5, 7.5e5,
  2.5e6, 7.5e6,
  2.5e7, 7.5e7
)
gather_breaks <- c(major_breaks_at, minor_breaks_at) %>% sort()

p <-
  ggplot(tb_p) +
  geom_abline(slope = 1, intercept = 0, color = alpha("#aaaaaa", 0.3)) +
  geom_abline(slope = 1, intercept = 2, color = alpha("#d7191c", 0.3)) +
  geom_abline(slope = 1, intercept = -2, color = alpha("#2c7bb6", 0.3)) +
  geom_point(
    aes(abun1, abun2, color = color, size = color),
    shape = 19, stroke = 0
  ) +
  # scale_x_continuous(
  #   trans = "log2",
  #   breaks = major_breaks_at,
  #   minor_breaks = minor_breaks_at
  # ) +
  # scale_y_continuous(
  #   trans = "log2",
  #   breaks = major_breaks_at,
  #   minor_breaks = minor_breaks_at
  # ) +
  scale_x_continuous(
    trans = "log2",
    breaks = major_breaks_at,
    minor_breaks = minor_breaks_at,
    labels = scales::number
  ) +
  scale_y_continuous(
    trans = "log2",
    breaks = major_breaks_at,
    minor_breaks = minor_breaks_at
  ) +
  scale_color_manual(
    values = c(
      "up" = "#d7191c",
      "not sign" = alpha("#aaaaaa", 1),
      "down" = "#2c7bb6"
    ),
    guide = guide_legend(
      title = "Ratio(JS005-HCP/JS005-FS)",
      override.aes = list(size = 2)
    )
  ) +
  scale_size_manual(
    values = c("up" = 2, "not sign" = 1.5, "down" = 2),
    guide = "none"
  ) +
  theme_bw() +
  theme(
    # panel.grid.minor = element_blank(),
    axis.text.x = element_text(hjust = 1, angle = 45)
  ) +
  coord_fixed() +
  labs(
    x = "JS005-FS abuandance (log2 scaled)",
    y = "JS005-HCP abuandance (log2 scaled)"
  )
p

ggsave(
  p, filename = path(a$pdo, "abun_and_ratio.pdf"),
  height = 6, width = 8
)
