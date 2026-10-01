suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
d <- "results"; out <- "figures"
FAM <- "Liberation Serif"; PT <- 8; TEAL <- "#4E7C8A"; ORANGE <- "#C97B3A"
base_theme <- theme_minimal(base_size = PT, base_family = FAM) +
  theme(text = element_text(size = PT, colour = "black"), axis.text = element_text(size = PT, colour = "black"),
        plot.title = element_text(size = PT, face = "bold", margin = margin(b = 2)),
        legend.text = element_text(size = PT), legend.title = element_text(size = PT),
        panel.grid.minor = element_blank(), plot.margin = margin(2, 6, 2, 2))

loo <- read_csv(file.path(d, "leave_one_out_influence.csv"), show_col_types = FALSE, col_types = cols(dropped = "c"))
full <- tibble(reference = 39.2, artifact = 16.0)
pA <- ggplot(loo, aes(reference, artifact)) +
  geom_point(data = full, shape = 4, size = 2.2, stroke = .8, colour = "black") +
  geom_point(colour = ifelse(loo$dropped == "47", ORANGE, TEAL), size = ifelse(loo$dropped == "47", 2.6, 1.5), alpha = .9) +
  annotate("text", x = 9.0, y = 7.5, label = "without participant 47", family = FAM, size = PT / .pt, hjust = 0, colour = ORANGE) +
  annotate("text", x = 39.2, y = 11.6, label = "all 71 (\u00d7)", family = FAM, size = PT / .pt, hjust = 0.5) +
  annotate("text", x = 47.5, y = 21.2, label = "any other participant removed", family = FAM, size = PT / .pt, hjust = 1, colour = TEAL) +
  scale_x_continuous(limits = c(0, 48), breaks = seq(0, 40, 10)) + scale_y_continuous(limits = c(0, 22), breaks = seq(0, 20, 5)) +
  labs(x = "Variance explained by reference (%)", y = "Variance explained by artifact handling (%)", title = "A. Leave one participant out") +
  base_theme

rec <- read_csv(file.path(d, "influential_recording.csv"), show_col_types = FALSE, col_types = cols(subject = "c")) %>%
  mutate(artifact = recode(artifact, ica = "Threshold + ICA", threshold = "Threshold only"),
         condition = factor(recode(condition, NS = "Normal\nsleep", SD = "Sleep\ndeprived"), levels = c("Normal\nsleep", "Sleep\ndeprived")))
band <- rec %>% filter(artifact == "Threshold only")
dodge <- position_dodge(width = 0.3)
pB <- ggplot(rec, aes(condition, theta_power, colour = artifact, group = artifact)) +
  geom_crossbar(data = band, aes(x = condition, y = others_median, ymin = others_p5, ymax = others_p95), inherit.aes = FALSE,
                width = 0.62, fill = "grey90", colour = "grey60", linewidth = .3, fatten = 1.2) +
  geom_line(linewidth = .7, position = dodge) +
  geom_errorbar(aes(ymin = theta_power - theta_sd, ymax = theta_power + theta_sd), width = .14, linewidth = .45, position = dodge, show.legend = FALSE) +
  geom_point(size = 2.2, position = dodge) +
  geom_text(data = rec %>% filter(condition == "Normal\nsleep"), aes(label = sprintf("%.1f", theta_power), y = theta_power),
            family = FAM, size = PT / .pt, hjust = 1.6, position = dodge, show.legend = FALSE) +
  annotate("text", x = 1.5, y = 1.38, label = "grey: other 70 participants\n(median, 5th\u201395th percentile)", family = FAM, size = 6 / .pt, hjust = 0.5, vjust = 0.5, colour = "grey35", lineheight = .8) +
  scale_colour_manual(values = c("Threshold + ICA" = TEAL, "Threshold only" = ORANGE), name = NULL) +
  scale_y_log10(limits = c(1, 120), breaks = c(1, 3, 10, 30, 100)) + scale_x_discrete(expand = expansion(add = 0.6)) +
  labs(x = NULL, y = expression(paste("Theta power (", mu, "V"^2, ", log scale)")), title = "B. Participant 47, centro-temporal") +
  base_theme + theme(legend.position = "bottom", legend.direction = "vertical", legend.key.size = unit(8, "pt"), legend.margin = margin(0, 0, 0, 0), panel.grid.major.x = element_blank())

six <- read_csv(file.path(d, "six_pipelines_own_vs_common.csv"), show_col_types = FALSE) %>%
  mutate(pipe = paste0(recode(reference, average = "Average", cz = "Cz", mastoid = "Mastoid"), " / ", recode(artifact, ica = "ICA", threshold = "threshold")),
         pipe = factor(pipe, levels = rev(pipe)), lab = paste0("n = ", n_own))
sixl <- six %>% pivot_longer(c(own_median_dz, common_median_dz), names_to = "sample", values_to = "dz") %>%
  mutate(sample = factor(sample, levels = c("own_median_dz", "common_median_dz"), labels = c("Pipeline's own sample", "Same 57 participants")))
pC <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_segment(data = six, aes(x = common_median_dz, xend = own_median_dz, y = pipe, yend = pipe), colour = "grey70", linewidth = .6) +
  geom_point(data = sixl, aes(dz, pipe, colour = sample), size = 2.4) +
  geom_text(data = six, aes(x = -0.1, y = pipe, label = lab), family = FAM, size = 6.5 / .pt, hjust = 0, colour = "grey30") +
  scale_colour_manual(values = c("Pipeline's own sample" = ORANGE, "Same 57 participants" = TEAL), name = NULL) +
  scale_x_continuous(limits = c(-0.1, 0.5), breaks = seq(0, 0.5, 0.1)) +
  labs(x = expression(paste("Median theta effect (Cohen's ", italic(d)[z], ")")), y = NULL, title = "C. Six 2-s pipelines, two samples") +
  base_theme + theme(legend.position = "bottom", legend.direction = "vertical", legend.key.size = unit(8, "pt"), legend.margin = margin(0, 0, 0, 0), panel.grid.major.y = element_blank())

fig <- pA + pB + pC + plot_layout(widths = c(1.15, 0.9, 1.15))
ggsave(file.path(out, "Fig1_sample_composition.png"), fig, width = 7.0, height = 2.9, dpi = 300, bg = "white", device = png, type = "cairo")
cat("done\n")
