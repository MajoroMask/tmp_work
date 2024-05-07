
setwd("~/work/tmp_work/20230407_forSYH")

library(rlang)
library(openxlsx)
library(ggplot2)

tb_p <-
  openxlsx::read.xlsx("./_i/20230407_105vs107_formatted.xlsx") %>%
  as_tibble()
p <-
  ggplot(tb_p) +
  geom_point(
    aes(x = pI, y = log10(MW), size = abundance, fill = group, alpha = group),
    shape = 21, stroke = 0
  ) +
  scale_fill_manual(values = c("#f56025", "#3366cc")) +
  scale_size_continuous(
    range = c(1, 25),
    guide = guide_legend(override.aes = list(alpha = 1, stroke = 1))
  ) +
  scale_alpha_manual(values = c("F105" = 0.8, "F107" = 0.385)) +
  theme_bw()
ggsave(p, filename = "./_o/20230407_bubble.pdf", width = 8, height = 4.5)
