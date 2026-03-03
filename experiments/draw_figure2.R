#!/usr/bin/env Rscript
# ============================================================
# draw_figure2.R — Generate Figure 2 (TPR/FPR bar plots)
# ============================================================
# 4-panel figure: Continuous (top) and Binary (bottom) traits,
# TPR (left) and FPR (right) for MAGMA, PASCAL, and LSEA.
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

# --- Constants ---
BARS_WIDTH <- 0.5
ERRORS_WIDTH <- 0.2

pastel_colors <- c("path_small" = "#b1decc", "path_medium" = "#a1cff0",
                   "path_big" = "#9598f0", "path_random" = "#ebbcdb")

path_labels <- c("path_small" = "Small (17 genes)",
                 "path_medium" = "Medium (69 genes)",
                 "path_big" = "Large (199 genes)",
                 "path_random" = "Random (FPR)")

model_levels <- c("MAGMA:\nlinreg", "MAGMA:\nmean", "MAGMA:\ntop", "PASCAL", "LSEA")

# --- Helper functions ---
get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

prepare_data <- function(filepath) {
  data <- read.csv(filepath)
  data$model <- recode(data$model,
                       "linreg" = "MAGMA:\nlinreg",
                       "mean" = "MAGMA:\nmean",
                       "top" = "MAGMA:\ntop")
  data$model <- factor(data$model, levels = model_levels)
  data$path <- factor(data$path, levels = names(pastel_colors))
  return(data)
}

make_tpr_plot <- function(data, title_text) {
  ggplot(data, aes(x = model, y = score, fill = path)) +
    geom_bar(stat = "identity", position = position_dodge(0.6), width = BARS_WIDTH, color = "black", linewidth = 0.3) +
    geom_point(position = position_dodge(0.6), size = 1) +
    geom_errorbar(aes(ymin = min, ymax = max), position = position_dodge(0.6), linewidth = 0.5, width = ERRORS_WIDTH) +
    scale_fill_manual(values = pastel_colors, labels = path_labels) +
    labs(title = title_text, y = "TPR", x = "") +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "none"
    )
}

make_fpr_plot <- function(data, title_text) {
  ggplot(data, aes(x = model, y = score, fill = path)) +
    geom_bar(stat = "identity", position = position_dodge(0.6), width = BARS_WIDTH / 3, color = "black", linewidth = 0.3) +
    geom_point(position = position_dodge(0.6), size = 1) +
    geom_errorbar(aes(ymin = min, ymax = max), position = position_dodge(0.6), linewidth = 0.5, width = ERRORS_WIDTH / 3) +
    scale_fill_manual(values = pastel_colors, labels = path_labels) +
    labs(title = title_text, y = "FPR", x = "") +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "none"
    )
}

# --- Read data ---
cont_tpr <- prepare_data(file.path(data_dir, "TPR_to_draw_LSEA.csv"))
cont_fpr <- prepare_data(file.path(data_dir, "FPR_to_draw_LSEA.csv"))
bin_tpr  <- prepare_data(file.path(data_dir, "binTPR_to_draw_LSEA.csv"))
bin_fpr  <- prepare_data(file.path(data_dir, "binFPR_to_draw_LSEA.csv"))

# --- Create 4 panels ---
p1 <- make_tpr_plot(cont_tpr, "Continuous: TPR")
p2 <- make_fpr_plot(cont_fpr, "Continuous: FPR")
p3 <- make_tpr_plot(bin_tpr, "Binary: TPR")
p4 <- make_fpr_plot(bin_fpr, "Binary: FPR")

# --- Shared legend (all 4 path types) ---
dummy_data <- data.frame(
  model = factor(rep("LSEA", 4), levels = model_levels),
  path = factor(names(pastel_colors), levels = names(pastel_colors)),
  score = rep(0.5, 4), min = rep(0.4, 4), max = rep(0.6, 4)
)
legend_plot <- ggplot(dummy_data, aes(x = model, y = score, fill = path)) +
  geom_bar(stat = "identity", position = position_dodge(0.6)) +
  scale_fill_manual(values = pastel_colors, labels = path_labels) +
  theme(legend.position = "bottom", legend.title = element_blank())
shared_legend <- get_only_legend(legend_plot)

# --- Combine: 2x2 grid + shared legend ---
p_grid <- grid.arrange(p1, p2, p3, p4, ncol = 2)
p_final <- grid.arrange(p_grid, shared_legend, nrow = 2, heights = c(10, 1))

# --- Save ---
pdf_out <- file.path(output_dir, "Figure2_TPR_FPR.pdf")
png_out <- file.path(output_dir, "Figure2_TPR_FPR.png")

ggsave(filename = pdf_out, plot = p_final, width = 10, height = 10, device = "pdf")
ggsave(filename = png_out, plot = p_final, width = 10, height = 10, dpi = 300, device = "png")

cat("Saved:", pdf_out, "\n")
cat("Saved:", png_out, "\n")
