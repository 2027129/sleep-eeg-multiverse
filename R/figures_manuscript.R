suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
res <- "results"
out <- "figures"
FAM <- "Liberation Serif"; PT <- 8
TEAL <- "#4E7C8A"; ORANGE <- "#C97B3A"; PLUM <- "#7B5EA7"; GREY <- "grey60"
REG <- c("Centro-temporal" = TEAL, "Parieto-occipital" = PLUM, "Frontal" = ORANGE)
base_theme <- theme_minimal(base_size = PT, base_family = FAM) +
  theme(text = element_text(size = PT, colour = "black"), axis.text = element_text(size = PT, colour = "black"),
        plot.title = element_text(size = PT, face = "bold", margin = margin(b = 2)),
        legend.text = element_text(size = PT), legend.title = element_text(size = PT),
        panel.grid.minor = element_blank(), plot.margin = margin(2, 4, 2, 2))
save_png <- function(p, stem, w, h) ggsave(file.path(out, paste0(stem, ".png")), p, width = w, height = h, dpi = 300, bg = "white", device = png, type = "cairo")
recode_roi <- function(x) recode(x, centro_temporal = "Centro-temporal", parieto_occipital = "Parieto-occipital", frontal = "Frontal")

B <- read_csv(file.path(res, "coreB_epoch2_specification_curve.csv"), show_col_types = FALSE)
spec_curve <- function(cv, title, labels_on = TRUE) {
  d <- cv %>% arrange(dz) %>% mutate(rank = row_number(), ci = qt(.975, n - 1) / sqrt(n),
                                     region = factor(recode_roi(roi), levels = names(REG)),
                                     sig = factor(ifelse(sig, "Significant after FDR correction (filled)", "Not significant (open)"),
                                                  levels = c("Significant after FDR correction (filled)", "Not significant (open)")))
  top <- ggplot(d, aes(rank, dz)) +
    geom_linerange(aes(ymin = dz - ci, ymax = dz + ci), colour = "grey85", linewidth = .3) +
    geom_point(aes(colour = region, shape = sig, fill = region), size = 1.3, stroke = .5) +
    scale_colour_manual(values = REG, name = "Region") +
    scale_fill_manual(values = REG, name = "Region") +
    scale_shape_manual(values = c(21, 1), name = NULL, drop = FALSE) +
    guides(fill = "none", colour = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 21, size = 2, fill = unname(REG))),
           shape = guide_legend(order = 2, nrow = 1, override.aes = list(size = 2, colour = "black", fill = c("black", "white")))) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
    geom_hline(yintercept = median(d$dz), linetype = "dotted", colour = "black", linewidth = .4) +
    coord_cartesian(ylim = c(-0.4, 1.15)) +
    labs(x = NULL, y = if (labels_on) expression(paste("Effect size (Cohen's ", italic(d)[z], ")")) else NULL, title = title) +
    base_theme + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  bottom <- d %>% select(rank, roi, normalise, reference) %>%
    pivot_longer(-rank, names_to = "node", values_to = "option") %>%
    mutate(node = factor(node, levels = c("roi", "normalise", "reference"), labels = c("Region", "Normalization", "Reference")),
           option = recode(option, rel_4_30 = "Relative to 4–30 Hz", rel_1_30 = "Relative to 1–30 Hz", absolute = "Absolute",
                           centro_temporal = "Centro-temporal", parieto_occipital = "Parieto-occipital", frontal = "Frontal",
                           average = "Average", cz = "Cz", mastoid = "Mastoid")) %>%
    ggplot(aes(rank, option)) + geom_point(shape = "|", size = 1.6) +
    facet_wrap(~ node, ncol = 1, scales = "free_y") +
    labs(x = "Specifications ranked by effect size", y = NULL) +
    base_theme + theme(strip.text = element_text(hjust = 0, face = "italic", size = PT, margin = margin(1, 0, 1, 0)),
                       panel.grid = element_blank(), axis.text.x = element_blank(),
                       axis.text.y = if (labels_on) element_text(size = PT, colour = "black") else element_blank(),
                       panel.spacing.y = unit(1, "pt"))
  list(top = top, bottom = bottom)
}
pA <- spec_curve(B %>% filter(band == "theta", outcome == "mean"), "A. Theta mean power")
pB <- spec_curve(B %>% filter(band == "theta", outcome == "sd"), "B. Theta epoch-to-epoch variability", labels_on = FALSE)
fig1 <- wrap_plots(pA$top, pB$top, pA$bottom, pB$bottom, ncol = 2, heights = c(2.0, 1.8)) +
  plot_layout(guides = "collect") & theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0), legend.box = "vertical",
                                          legend.spacing.y = unit(0, "pt"), legend.justification = "left")
save_png(fig1, "Fig2_specification_curves", 7.0, 4.6)

vd <- read_csv(file.path(res, "variance_decomposition_all_sets.csv"), show_col_types = FALSE)
vd_plot <- vd %>% filter(set %in% c("unrestricted_theta_mean", "coreB_theta_mean"), term != "Residuals") %>%
  mutate(set = recode(set, unrestricted_theta_mean = "Own sample per pipeline (n = 14–71)",
                      coreB_theta_mean = "Same 57 participants"),
         term = recode(term, epoch_sec = "Epoch length", roi = "Region", artifact = "Artifact handling",
                       normalise = "Normalization", aperiodic = "Aperiodic removal", reference = "Reference"),
         term = factor(term, levels = c("Reference", "Artifact handling", "Region", "Normalization", "Aperiodic removal", "Epoch length")))
lv <- c("Own sample per pipeline (n = 14–71)", "Same 57 participants")
vd_plot <- vd_plot %>% bind_rows(tibble(term = factor("Epoch length", levels = levels(vd_plot$term)), set = lv[2], pct_variance = 0)) %>%
  mutate(set = factor(set, levels = rev(lv)), lab = ifelse(pct_variance > 0, sprintf("%.1f", pct_variance), ""))
fig2 <- ggplot(vd_plot, aes(fct_rev(term), pct_variance, fill = set)) +
  geom_col(position = position_dodge(width = .8), width = .75) +
  geom_text(aes(label = lab, group = set), family = FAM, size = PT / .pt, position = position_dodge(width = .8), hjust = -0.15) +
  scale_fill_manual(values = setNames(c(TEAL, ORANGE), rev(lv)), breaks = lv, name = NULL) +
  scale_y_continuous(limits = c(0, 46), expand = expansion(mult = c(0, .02))) + coord_flip() +
  labs(y = "Share of variance in theta mean effect size (%)", x = NULL, title = "A. Share of variance by decision node") +
  base_theme + theme(legend.position = "bottom", legend.direction = "vertical", legend.key.size = unit(8, "pt"), panel.grid.major.y = element_blank())
save_png(fig2, "Fig3A_variance", 3.35, 3.0)

cat("done\n")
