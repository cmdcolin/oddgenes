#!/usr/bin/env Rscript
# chry_composite.R — builds img/chry_composite.png
#
# Five panels on why an XX sample's chrY reads match a male sample's Y at 99.8%:
#   A  chrY ideogram with the regions that attract reads
#   B  chrY depth for NA12882 (XY) above his mother NA12878 (XX), with their
#      PAR1 reads realigned against a reference that keeps chrY PAR
#   C  identity scored as differing positions, on the 546,417 positions she covers
#   D  the coverage that percentage was computed on
#   E  MSY depth, the measurement that does separate XX from XY
#
# Samples are CEPH pedigree 1463 from the Platinum Pedigree Consortium, DRAGEN
# v4.2.4 against GRCh38: NA12878 (XX), her sons NA12882 and NA12886, her husband
# NA12877, his father NA12889 (the Y her sons carry), her father NA12891 (an
# unrelated Y), her daughter NA12879.
#
#   usage:  Rscript chry_composite.R [data_dir] [output.png]
#
# make_inputs.sh regenerates data_dir from the public BAMs; the counts quoted in
# panels C and D come from analysis.py and are hard-coded below.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
DATA <- if (length(args) >= 1) args[1] else "data"
OUT  <- if (length(args) >= 2) args[2] else "chry_composite.png"

CHRY <- 57227415 / 1e6
FLOOR <- 0.004; YTOP <- 120

SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRIDC <- "#e1e0d9"; AXISC <- "#c3c2b7"
MALE <- "#2a78d6"; FEMALE <- "#eb6834"
MALE_L <- "#a9c9f0"; FEMALE_L <- "#f6b79a"

BARW <- 3.8

base_theme <- theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = GRIDC, linewidth = 0.3),
        axis.text = element_text(colour = MUTED, size = 9.6),
        axis.title = element_text(colour = INK2, size = 10),
        plot.title = element_text(colour = INK, size = 12, face = "plain", hjust = 0),
        plot.subtitle = element_text(colour = MUTED, size = 9.4, hjust = 0),
        plot.margin = margin(3, 4, 3, 4),
        legend.position = "none")

## ---------------------------------------------------------------- where the reads land
norm <- read.delim(file.path(DATA, "counts.tsv"), header = FALSE,
                   col.names = c("sample", "reads", "q30", "auto"))
wdr <- read.delim(file.path(DATA, "windows.tsv"))

LV <- c("son", "son · PAR-kept", "mother", "mother · PAR-kept")

parbins <- function(f, sample) {
  a <- norm$auto[norm$sample == sample]
  read.delim(gzfile(file.path(DATA, f)), header = FALSE, col.names = c("chr","start","end","dep")) %>%
    filter(chr == "chrY", start >= 300000, start < 2781479) %>%
    transmute(mb = start / 1e6, rel = pmax(dep / a, FLOOR))
}

cov <- bind_rows(
  wdr %>% transmute(mb = start/1e6,
                    rel = ifelse(nfrac >= 0.5, NA_real_, pmax(NA12882_son, FLOOR)),
                    series = LV[1], panel = "XY  NA12882 son"),
  parbins("NA12882_son.parwin.regions.bed.gz", "NA12882_son") %>%
    mutate(series = LV[2], panel = "XY  NA12882 son"),
  wdr %>% transmute(mb = start/1e6,
                    rel = ifelse(nfrac >= 0.5, NA_real_, pmax(NA12878_mother, FLOOR)),
                    series = LV[3], panel = "XX  NA12878 mother"),
  parbins("NA12878_mother.parwin.regions.bed.gz", "NA12878_mother") %>%
    mutate(series = LV[4], panel = "XX  NA12878 mother")
) %>% mutate(series = factor(series, levels = LV),
             panel = factor(panel, levels = c("XY  NA12882 son", "XX  NA12878 mother")))

REGIONS <- tribble(
  ~start, ~end,   ~label,           ~kind,
  0.01,   2.781,  "PAR1",           "pseudoautosomal (masked here)",
  56.888, 57.217, "PAR2",           "pseudoautosomal (masked here)",
  3.050,  6.234,  "XTR",            "X-transposed (~99% like Xq21)",
  10.0,   12.5,   "beta satellite", "beta satellite array",
  56.673, 56.888, "",               "beta satellite array",
  26.673, 56.673, "Yq12",           "heterochromatin — 30 Mb of N in the reference")
KIND_COL <- c("pseudoautosomal (masked here)" = "#4a3aa7",
              "X-transposed (~99% like Xq21)" = "#1baf7a",
              "beta satellite array" = "#eda100",
              "heterochromatin — 30 Mb of N in the reference" = "#c3c2b7")

cyto <- read_tsv(file.path(DATA, "cytoband_chrY.tsv"), show_col_types = FALSE) %>%
  mutate(x1 = start / 1e6, x2 = end / 1e6)
STAIN <- c(gneg = "#f7f7f5", gpos25 = "#c9c8c2", gpos50 = "#9d9c95",
           gpos75 = "#6f6e68", gpos100 = "#3d3c38", acen = "#c0392b", gvar = "#dfe6ef")

pA <- ggplot() +
  geom_rect(data = cyto, aes(xmin = x1, xmax = x2, ymin = 0.55, ymax = 1.0, fill = stain),
            colour = AXISC, linewidth = 0.15) +
  scale_fill_manual(values = STAIN, guide = "none") +
  geom_rect(data = REGIONS, aes(xmin = start, xmax = end, ymin = 0.05, ymax = 0.40, colour = kind),
            fill = NA, linewidth = 0.9) +
  scale_colour_manual(values = KIND_COL, guide = "none") +
  geom_text(data = REGIONS %>% filter(label != "") %>%
              mutate(lx = ifelse(label == "PAR2", start - 0.4, (start + end) / 2),
                     lh = ifelse(label == "PAR2", 1, 0.5)),
            aes(x = lx, y = -0.25, label = label, hjust = lh), size = 3.2, colour = INK2) +
  scale_x_continuous(limits = c(0, CHRY), expand = c(0.005, 0)) +
  scale_y_continuous(limits = c(-0.5, 1.05), breaks = c(0, 0.5, 1), labels = c("0.01", "0.5", "100")) +
  labs(title = "A · chrY and the regions that attract reads",
       subtitle = "PAR masked here · XTR ~99% like Xq21 · Yq12 is 30 Mb of N",
       y = "Y depth / autosomal depth", x = NULL) +
  base_theme +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(colour = NA, size = 8.5),
        axis.title.y = element_text(colour = NA, size = 9), panel.grid = element_blank())

parlab <- data.frame(mb = 3.4, rel = 25, panel = factor("XY  NA12882 son", levels = levels(cov$panel)))

pB <- ggplot(cov, aes(x = mb)) +
  geom_ribbon(aes(ymin = FLOOR, ymax = rel, fill = series), outline.type = "upper") +
  geom_text(data = parlab, aes(x = mb, y = rel), inherit.aes = FALSE, hjust = 0, size = 3.2,
            colour = MUTED, label = "pale = PAR1 remapped to a reference that keeps PAR") +
  geom_hline(yintercept = 0.5, colour = AXISC, linetype = "22", linewidth = 0.4) +
  annotate("text", x = CHRY, y = 0.55, label = "one copy", hjust = 1, vjust = 0, size = 3.2, colour = MUTED) +
  facet_wrap(~panel, ncol = 1) +
  scale_fill_manual(values = setNames(c(MALE, MALE_L, FEMALE, FEMALE_L), LV), guide = "none") +
  scale_y_log10(limits = c(FLOOR, YTOP), breaks = c(0.01, 0.1, 0.5, 10, 100),
                labels = c("0.01", "0.1", "0.5", "10", "100")) +
  scale_x_continuous(limits = c(0, CHRY), expand = c(0.005, 0)) +
  labs(title = "B · chrY reads, son and mother",
       y = "Y depth / autosomal depth", x = "chrY position (Mb)") +
  base_theme +
  theme(strip.text = element_text(colour = INK, size = 10.5, hjust = 0))

## ---------------------------------------------------------------- what a % match measures
d <- data.frame(
  who = c("paternal grandfather & grandson · same Y",
          "maternal grandfather & grandson · different Ys",
          "mother & son · no Y in common",
          "mother & maternal grandfather · no Y in common"),
  diff = c(6, 72, 887, 893),
  pct  = c(99.9989, 99.9864, 99.8375, 99.8317),
  fill = c(MALE, MALE, FEMALE, FEMALE))
d$who <- factor(d$who, levels = rev(d$who))
pC <- ggplot(d, aes(diff, who, colour = fill)) +
  geom_segment(aes(x = 4.6, xend = diff, yend = who), linewidth = BARW, lineend = "butt") +
  geom_text(aes(label = sprintf("%s  ·  %s", format(diff, big.mark = ","), sprintf("%.4f%%", pct))),
            hjust = -0.15, size = 3.4, colour = INK2) +
  scale_colour_identity() +
  scale_x_log10(limits = c(4.4, 20000), breaks = c(10, 100, 1000),
                labels = c("10", "100", "1,000"), expand = c(0, 0)) +
  labs(title = "C · the same matches, counted as differences",
       subtitle = "over the 546,417 positions the mother covers",
       x = "differing positions (log scale)", y = NULL) +
  base_theme + theme(axis.text.y = element_text(size = 9, colour = INK2))

e <- data.frame(
  who = factor(c("XX  mother", "XY  son"), levels = c("XY  son", "XX  mother")),
  pct = c(7.8, 100), fill = c(FEMALE, MALE))
pD <- ggplot(e, aes(pct, who, colour = fill)) +
  geom_segment(aes(x = 0, xend = pct, yend = who), linewidth = BARW, lineend = "butt") +
  geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.15, size = 3.4, colour = INK2) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 138), breaks = c(0, 50, 100), expand = c(0, 0)) +
  labs(title = "D · the coverage it was computed on",
       subtitle = "% of the same 7 Mb with any read",
       x = NULL, y = NULL) +
  base_theme + theme(axis.text.y = element_text(size = 9, colour = INK2))

w <- read.delim(file.path(DATA, "windows.tsv")) %>% filter(nfrac < 0.5, class == "MSY unique")
g <- w %>%
  select(NA12878_mother, NA12879_daughter, NA12882_son, NA12886_son2, NA12877_father,
         NA12889_unrel_male, NA12891_her_father) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "depth") %>%
  mutate(karyo = ifelse(grepl("mother|daughter", sample), "XX", "XY"),
         sample = sub("_.*", "", sample), depth = pmax(depth, 0.001))
g$sample <- factor(g$sample, levels = rev(c("NA12878","NA12879","NA12882","NA12886","NA12877","NA12889","NA12891")))
pE <- ggplot(g, aes(depth, sample, colour = karyo)) +
  geom_jitter(height = 0.18, size = 0.25, alpha = 0.25) +
  stat_summary(fun = median, geom = "point", size = 2.6, shape = 18) +
  geom_vline(xintercept = 0.5, linetype = "22", colour = AXISC, linewidth = 0.4) +
  annotate("text", x = Inf, y = 7.8, label = "one Y present", hjust = 1, size = 3.4, colour = INK2) +
  annotate("text", x = 0.00105, y = 7.8, label = "mapping artifact only", hjust = 0, size = 3.4, colour = INK2) +
  scale_colour_manual(values = c(XX = FEMALE, XY = MALE)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.5, 1), labels = c("0.001","0.01","0.1","0.5","1"),
                expand = expansion(mult = c(0.02, 0.06))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.25)),
                   labels = c(NA12878 = "XX  NA12878 mother", NA12879 = "XX  NA12879 daughter",
                              NA12882 = "XY  NA12882 son", NA12886 = "XY  NA12886 son",
                              NA12877 = "XY  NA12877 father", NA12889 = "XY  NA12889 paternal grandfather",
                              NA12891 = "XY  NA12891 maternal grandfather")) +
  coord_cartesian(clip = "off") +
  labs(title = "E · what does separate them: Y-specific depth",
       subtitle = "10 kb windows with no chrX counterpart · diamond = median",
       x = "Y depth / autosomal depth", y = NULL) + base_theme

left  <- pA / pB + plot_layout(heights = c(0.58, 2.55))
right <- pC / pD / pE + plot_layout(heights = c(1.2, 0.62, 2.3))


## ---------------------------------------------------------------- the logical error
errtxt <- data.frame(
  y = c(4, 3, 2, 1),
  n = c("1", "2", "3", "4"),
  t = c("Computed only where the mother has reads — 7.8% of the Y. The 92% with no coverage IS the evidence, and the statistic discards it.",
        "No control. An unrelated man's Y scores 99.832% against the mother, the son's 99.837%: the positions both cover are X-homologous sequence, where any two people match.",
        "99.8% is not a high score. The paternal grandfather and grandson share a Y and differ at just 6 positions — errors or de novo mutations. Mother and son differ at 887.",
        "Depth is the measurement with power: a Y/autosomal depth ratio of 0.004 means no Y, 0.5 means one copy, and nothing lands in the hundred-fold gap between."))

pErr <- ggplot(errtxt) +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0.4, ymax = 4.6, fill = "#f4efe6", colour = "#e0d6c2", linewidth = 0.4) +
  geom_text(aes(x = 0.16, y = y, label = n), hjust = 0, size = 4.2, colour = "#b5651d", fontface = "bold") +
  geom_text(aes(x = 0.42, y = y, label = t), hjust = 0, size = 3.9, colour = INK) +
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 5.35), expand = c(0, 0)) +
  labs(title = "The four mistakes") +
  theme_void(base_size = 12) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        plot.title = element_text(colour = INK, size = 12.5, hjust = 0, margin = margin(0, 0, 4, 4)),
        plot.margin = margin(4, 4, 2, 4))

top <- (wrap_elements(full = left) | wrap_elements(full = right)) +
  plot_layout(widths = c(1, 1.18))
fig <- (top / pErr) + plot_layout(heights = c(1, 0.26)) +
  plot_annotation(
    title = "chrY reads in XX samples: where they come from, and why \"99.9% identical\" measures nothing",
    subtitle = "CEPH pedigree 1463 · Platinum Pedigree DRAGEN v4.2.4 on GRCh38 · depth ratios from 10 kb windows over a 5 Mb chr20 baseline",
    theme = theme(plot.title = element_text(size = 12, colour = INK),
                  plot.subtitle = element_text(size = 8, colour = MUTED),
                  plot.background = element_rect(fill = SURFACE, colour = NA)))

ggsave(OUT, fig, width = 14.4, height = 9.0, dpi = 172, bg = SURFACE)
cat("wrote chry_composite.png\n")
