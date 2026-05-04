#!/usr/bin/env Rscript
# ============================================================
# draw_figure2.R — Generate Figure 2 (TPR/FPR bar plots)
# ============================================================
# 2x2 panel figure: rows = Continuous / Binary,
#                    cols = TPR (3 pathway sizes) / FPR (random).
# X-axis  = pathway size, fill color = tool.
# Right column (FPR) is narrower than left (TPR).
#
# Usage:
#   Rscript draw_figure2.R <data_dir> <output_dir>
#
# Expected data files in <data_dir>:
#   TPR_to_draw_LSEA.csv, FPR_to_draw_LSEA.csv
#   binTPR_to_draw_LSEA.csv, binFPR_to_draw_LSEA.csv
# ============================================================

library(ggplot2)
library(dplyr)
library(gridExtra)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript draw_figure2.R <data_dir> <output_dir>\n")
  quit(status = 1)
}
data_dir <- args[1]
output_dir <- args[2]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Tool color palette (MAGMA = purple family, PASCAL/LSEA = distinct) ---
tool_colors <- c("MAGMA: linreg" = "#9E82C0",
                 "MAGMA: mean"   = "#C8B0E0",
                 "MAGMA: top"    = "#7558A0",
                 "PASCAL"        = "#5BC0A0",
                 "LSEA"          = "#E8887A")

path_labels <- c("path_small"  = "Small\n(17 genes)",
                 "path_medium" = "Medium\n(69 genes)",
                 "path_big"    = "Large\n(199 genes)",
                 "path_random" = "Random\n ")

# --- Helper functions ---
get_only_legend <- function(plot) {
  plot_table  <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  plot_table$grobs[[legend_plot]]
}

prepare_data <- function(filepath) {
  data <- read.csv(filepath)
  data$model <- recode(data$model,
                       "linreg" = "MAGMA: linreg",
                       "mean"   = "MAGMA: mean",
                       "top"    = "MAGMA: top")
  data$model <- factor(data$model, levels = names(tool_colors))
  data$path  <- factor(data$path,
                       levels = c("path_small", "path_medium",
                                  "path_big", "path_random"))
  return(data)
}

shared_theme <- theme_minimal() +
  theme(
    legend.title       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    axis.text.x        = element_text(size = 10),
    axis.text.y        = element_text(size = 11),
    axis.title.y       = element_text(size = 13),
    plot.title         = element_text(size = 14, face = "bold"),
    legend.position    = "none",
    plot.margin        = margin(t = 5, r = 5, b = 30, l = 5, unit = "pt")
  )

make_plot <- function(data, title_text, ylab_text) {
  p <- ggplot(data, aes(x = path, y = score, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(0.7),
             width = 0.6, color = "#444444", linewidth = 0.3) +
    geom_errorbar(aes(ymin = min, ymax = max),
                  position = position_dodge(0.7),
                  linewidth = 0.4, width = 0.15) +
    scale_fill_manual(values = tool_colors) +
    scale_x_discrete(labels = path_labels) +
    labs(title = title_text, y = ylab_text, x = "") +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2),
                       expand = c(0, 0)) +
    shared_theme
  return(p)
}

# --- Read data ---
cont_tpr <- prepare_data(file.path(data_dir, "TPR_to_draw_LSEA.csv"))
cont_fpr <- prepare_data(file.path(data_dir, "FPR_to_draw_LSEA.csv"))
bin_tpr  <- prepare_data(file.path(data_dir, "binTPR_to_draw_LSEA.csv"))
bin_fpr  <- prepare_data(file.path(data_dir, "binFPR_to_draw_LSEA.csv"))

# --- Create 4 panels (x = pathway size, fill = tool) ---
p1 <- make_plot(cont_tpr, "Continuous: TPR", "TPR")
p2 <- make_plot(cont_fpr, "Continuous: FPR", "FPR")
p3 <- make_plot(bin_tpr,  "Binary: TPR", "TPR")
p4 <- make_plot(bin_fpr,  "Binary: FPR", "FPR")

# --- Shared legend (horizontal, at bottom) ---
dummy_data <- data.frame(
  model = factor(names(tool_colors), levels = names(tool_colors)),
  path  = factor(rep("path_small", 5),
                 levels = c("path_small", "path_medium",
                            "path_big", "path_random")),
  score = rep(0.5, 5), min = rep(0.4, 5), max = rep(0.6, 5)
)
legend_plot <- ggplot(dummy_data, aes(x = path, y = score, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge(0.6),
           color = "#444444", linewidth = 0.4) +
  scale_fill_manual(values = tool_colors) +
  guides(fill = guide_legend(override.aes = list(color = "#444444",
                                                  linewidth = 0.5))) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(0.7, "cm"),
    legend.spacing.x   = unit(0.3, "cm")
  )
shared_legend <- get_only_legend(legend_plot)

# --- Combine: 2x2 grid (right column narrower) + shared legend ---
p_grid  <- grid.arrange(p1, p2, p3, p4,
                         layout_matrix = rbind(c(1, 2), c(3, 4)),
                         widths = c(3, 1.3))
p_final <- grid.arrange(p_grid, shared_legend, nrow = 2, heights = c(10, 1.2))

# --- Save ---
pdf_out <- file.path(output_dir, "Figure2_TPR_FPR.pdf")
png_out <- file.path(output_dir, "Figure2_TPR_FPR.png")

ggsave(filename = pdf_out, plot = p_final, width = 11, height = 10,
       device = "pdf", bg = "white")
ggsave(filename = png_out, plot = p_final, width = 11, height = 10,
       dpi = 300, device = "png", bg = "white")

cat("Saved:", pdf_out, "\n")
cat("Saved:", png_out, "\n")
