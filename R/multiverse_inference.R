suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
})

set.seed(20260927)
N_PERM   <- 2000
out_dir  <- "results"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

spec <- c("reference", "artifact", "epoch_sec", "normalise",
          "aperiodic", "outcome", "band", "roi")

feat <- read_parquet(file.path(out_dir, "features_all_pipelines.parquet"))
cat("Feature rows:", nrow(feat), " participants:", n_distinct(feat$subject), "\n")

retained <- feat %>%
  distinct(subject, condition, reference, artifact, epoch_sec) %>%
  count(reference, artifact, epoch_sec, name = "recordings") %>%
  pivot_wider(names_from = epoch_sec, values_from = recordings,
              names_prefix = "epoch_")
write_csv(retained, file.path(out_dir, "sample_retained_table.csv"))
print(retained)

presence <- feat %>% distinct(subject, condition, reference, artifact, epoch_sec)

complete_B <- presence %>% filter(epoch_sec == 2) %>%
  count(subject) %>% filter(n == 3 * 2 * 2) %>% pull(subject)
complete_A <- presence %>% filter(epoch_sec %in% c(2, 4)) %>%
  count(subject) %>% filter(n == 3 * 2 * 2 * 2) %>% pull(subject)

cat("Core B (2 s) participants:", length(complete_B),
    "   Core A (2 & 4 s) participants:", length(complete_A), "\n")

core_B <- feat %>% filter(epoch_sec == 2, subject %in% complete_B)
core_A <- feat %>% filter(epoch_sec %in% c(2, 4), subject %in% complete_A)

diff_matrix <- function(core) {
  w <- core %>%
    select(all_of(spec), subject, condition, value) %>%
    pivot_wider(names_from = condition, values_from = value) %>%
    mutate(diff = SD - NS)
  specs <- w %>% distinct(across(all_of(spec))) %>%
    arrange(across(all_of(spec))) %>% mutate(spec_id = row_number())
  w <- w %>% inner_join(specs, by = spec)
  M <- w %>% select(spec_id, subject, diff) %>%
    pivot_wider(names_from = subject, values_from = diff) %>% arrange(spec_id)
  mat <- as.matrix(M[, -1]); rownames(mat) <- M$spec_id
  bad <- !is.finite(mat)
  if (any(bad)) {
    warning(sum(bad), " non-finite differences; affected specifications dropped")
    keep  <- rowSums(bad) == 0
    mat   <- mat[keep, , drop = FALSE]
    specs <- specs %>% filter(spec_id %in% as.integer(rownames(mat)))
  }
  list(specs = specs, mat = mat)
}

curve_from_matrix <- function(mat) {
  n  <- ncol(mat)
  m  <- rowMeans(mat)
  s  <- sqrt(rowSums(sweep(mat, 1, m)^2) / (n - 1))
  dz <- m / s
  tt <- dz * sqrt(n)
  p  <- 2 * pt(abs(tt), df = n - 1, lower.tail = FALSE)
  tibble(spec_id = as.integer(rownames(mat)), n = n, dz = dz, p = p)
}

build_curve <- function(core, label) {
  dm <- diff_matrix(core)
  cv <- curve_from_matrix(dm$mat) %>%
    inner_join(dm$specs, by = "spec_id") %>%
    mutate(p_fdr = p.adjust(p, "BH"), sig = p_fdr < 0.05) %>%
    select(all_of(spec), n, dz, p, p_fdr, sig)
  write_csv(cv, file.path(out_dir, paste0(label, "_specification_curve.csv")))
  list(curve = cv, dm = dm)
}

B <- build_curve(core_B, "coreB_epoch2")
A <- build_curve(core_A, "coreA_epoch2and4")

summarise_curve <- function(cv, by = c("band", "outcome")) {
  cv %>% group_by(across(all_of(by))) %>%
    summarise(N = n(), median_dz = median(dz),
              q25 = quantile(dz, .25), q75 = quantile(dz, .75),
              min = min(dz), max = max(dz),
              pct_sig_fdr = 100 * mean(sig), pct_positive = 100 * mean(dz > 0),
              .groups = "drop")
}

cat("\nCORE B (primary): by band x outcome\n")
print(summarise_curve(B$curve), n = 20)
theta_mean_B <- B$curve %>% filter(band == "theta", outcome == "mean")
for (node in c("roi", "normalise", "reference", "artifact", "aperiodic")) {
  cat("\ntheta mean by", node, "\n"); print(summarise_curve(theta_mean_B, node))
}
cat("\nCORE A (sensitivity): by band x outcome\n")
print(summarise_curve(A$curve), n = 20)

decompose <- function(cv, formula) {
  a <- anova(lm(formula, data = cv %>% mutate(epoch_sec = factor(epoch_sec))))
  tibble(term = rownames(a), df = a$Df, sumsq = a$`Sum Sq`) %>%
    mutate(pct_variance = 100 * sumsq / sum(sumsq)) %>%
    arrange(desc(pct_variance))
}

vd_all_B   <- decompose(B$curve,
               dz ~ reference + artifact + normalise + aperiodic + outcome + band + roi)
vd_int_B   <- decompose(B$curve,
               dz ~ (reference + artifact + normalise + aperiodic + outcome + band + roi)^2)
vd_theta_B <- decompose(theta_mean_B, dz ~ reference + artifact + normalise + aperiodic + roi)
vd_theta_A <- decompose(A$curve %>% filter(band == "theta", outcome == "mean"),
               dz ~ reference + artifact + epoch_sec + normalise + aperiodic + roi)

unrestricted <- read_csv(file.path(out_dir, "specification_curve.csv"), show_col_types = FALSE)
vd_theta_U <- decompose(unrestricted %>% filter(band == "theta", outcome == "mean"),
               dz ~ reference + artifact + epoch_sec + normalise + aperiodic + roi)

cat("\nVariance decomposition, core B, all specifications (main effects)\n"); print(vd_all_B)
cat("\nCore B, two-way interactions (top 10)\n"); print(head(vd_int_B, 10))
cat("\nCore B, theta-mean subset\n"); print(vd_theta_B)
cat("\nUnrestricted multiverse, theta-mean subset (for contrast)\n"); print(vd_theta_U)

bind_rows(vd_all_B %>% mutate(set = "coreB_all"),
          vd_int_B %>% mutate(set = "coreB_all_2way"),
          vd_theta_B %>% mutate(set = "coreB_theta_mean"),
          vd_theta_A %>% mutate(set = "coreA_theta_mean"),
          vd_theta_U %>% mutate(set = "unrestricted_theta_mean")) %>%
  write_csv(file.path(out_dir, "variance_decomposition_all_sets.csv"))

joint_tests <- function(mat, expected_sign = 1, B = N_PERM) {
  n <- ncol(mat)
  stat <- function(cv) {
    z <- qnorm(1 - cv$p / 2) * sign(cv$dz)
    c(median_dz          = expected_sign * median(cv$dz),
      share_sig_expected = mean(cv$p < 0.05 & sign(cv$dz) == expected_sign),
      mean_z             = expected_sign * mean(z))
  }
  obs  <- stat(curve_from_matrix(mat))
  null <- t(replicate(B, {
    flips <- sample(c(-1, 1), n, replace = TRUE)
    stat(curve_from_matrix(sweep(mat, 2, flips, `*`)))
  }))
  p <- sapply(names(obs), function(k) (1 + sum(null[, k] >= obs[[k]])) / (B + 1))
  tibble(statistic   = names(obs),
         observed    = unname(obs),
         null_median = apply(null, 2, median),
         null_q95    = apply(null, 2, quantile, .95),
         p_perm      = unname(p))
}

families <- tribble(
  ~band,   ~outcome, ~expected_sign,
  "theta", "mean",    1,
  "theta", "sd",      1,
  "theta", "cv",      1,
  "alpha", "mean",   -1,
  "alpha", "sd",      1,
  "alpha", "cv",      1,
  "beta",  "mean",    1,
  "beta",  "sd",      1,
  "beta",  "cv",      1
)

cat("\nRunning", N_PERM, "permutations per family ...\n")
perm_results <- pmap_dfr(families, function(band, outcome, expected_sign) {
  ids <- B$dm$specs %>%
    filter(band == !!band, outcome == !!outcome) %>% pull(spec_id)
  joint_tests(B$dm$mat[as.character(ids), , drop = FALSE], expected_sign) %>%
    mutate(band = band, outcome = outcome, expected_sign = expected_sign,
           n_specs = length(ids), .before = 1)
})
write_csv(perm_results, file.path(out_dir, "permutation_joint_tests_coreB.csv"))
cat("\nPermutation joint tests, core B\n"); print(perm_results, n = 40)

TEAL <- "#4E7C8A"

spec_curve_plot <- function(cv, title, ylab_on = TRUE) {
  d <- cv %>% arrange(dz) %>%
    mutate(rank = row_number(), ci = qt(.975, n - 1) / sqrt(n))
  top <- ggplot(d, aes(rank, dz)) +
    geom_linerange(aes(ymin = dz - ci, ymax = dz + ci), colour = "grey80", linewidth = .4) +
    geom_point(aes(colour = sig), size = 1.1) +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = TEAL),
                        labels = c("n.s.", "p_FDR < 0.05"), name = NULL) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = median(d$dz), linetype = "dotted", colour = TEAL) +
    coord_cartesian(ylim = c(-0.4, 1.15)) +
    labs(x = NULL, y = if (ylab_on) "SD \u2212 NS (Cohen's dz)" else NULL, title = title,
         subtitle = sprintf("%d specifications, n = %d each \u00b7 median dz = %.2f\n%.0f%% significant after FDR \u00b7 %.0f%% positive",
                            nrow(d), d$n[1], median(d$dz), 100 * mean(d$sig), 100 * mean(d$dz > 0))) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(), legend.position = "top",
          panel.grid.minor = element_blank())
  bottom <- d %>% select(rank, roi, normalise, reference) %>%
    pivot_longer(-rank, names_to = "node", values_to = "option") %>%
    mutate(node   = factor(node, levels = c("roi", "normalise", "reference"),
                           labels = c("region", "normalisation", "reference")),
           option = recode(option, rel_4_30 = "relative 4\u201330 Hz", rel_1_30 = "relative 1\u201330 Hz",
                           centro_temporal = "centro-temporal", parieto_occipital = "parieto-occipital",
                           cz = "Cz")) %>%
    ggplot(aes(rank, option)) + geom_point(shape = "|", size = 2) +
    facet_grid(node ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "Specifications, sorted by effect size", y = NULL) +
    theme_minimal(base_size = 8) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0),
          panel.grid = element_blank())
  list(top = top, bottom = bottom)
}

p_mean <- spec_curve_plot(theta_mean_B, "A. Theta mean power")
p_sd   <- spec_curve_plot(B$curve %>% filter(band == "theta", outcome == "sd"),
                          "B. Theta epoch-to-epoch variability", ylab_on = FALSE)

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig1 <- ((p_mean$top / p_mean$bottom) + plot_layout(heights = c(2.2, 1.4))) |
          ((p_sd$top / p_sd$bottom) + plot_layout(heights = c(2.2, 1.4)))
  ggsave(file.path(fig_dir, "Fig1_specification_curves_coreB.png"), fig1,
         width = 11, height = 6.8, dpi = 300, bg = "white")
} else {
  ggsave(file.path(fig_dir, "Fig1A_theta_mean_curve.png"),  p_mean$top, width = 6, height = 4, dpi = 200, bg = "white")
  ggsave(file.path(fig_dir, "Fig1A_theta_mean_nodes.png"),  p_mean$bottom, width = 6, height = 2.5, dpi = 200, bg = "white")
  ggsave(file.path(fig_dir, "Fig1B_theta_sd_curve.png"),    p_sd$top, width = 6, height = 4, dpi = 200, bg = "white")
  ggsave(file.path(fig_dir, "Fig1B_theta_sd_nodes.png"),    p_sd$bottom, width = 6, height = 2.5, dpi = 200, bg = "white")
  message("patchwork not installed: Figure 1 saved as four separate panels")
}

vd_plot <- bind_rows(
    vd_theta_U %>% mutate(set = "All specifications: 432 theta-mean (n = 14\u201371)"),
    vd_theta_B %>% mutate(set = "Fixed sample: 108 theta-mean (2-s epochs, n = 57)")) %>%
  filter(term != "Residuals") %>%
  mutate(term = recode(term, epoch_sec = "epoch length", roi = "region (ROI)",
                       artifact = "artifact handling", normalise = "normalisation",
                       aperiodic = "aperiodic removal"),
         term = factor(term, levels = c("reference", "artifact handling", "region (ROI)",
                                        "normalisation", "aperiodic removal", "epoch length")))

fig2 <- ggplot(vd_plot, aes(pct_variance, fct_rev(term), fill = set)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_text(aes(label = sprintf("%.1f", pct_variance)),
            position = position_dodge(width = .8), hjust = -0.2, size = 2.6) +
  scale_fill_manual(values = c("grey75", TEAL), name = NULL) +
  labs(x = "% of variance in theta-mean effect size (main effects)", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "Fig2_variance_decomposition.png"), fig2, width = 8, height = 4, dpi = 300, bg = "white")

cat("\nDone. Outputs in", out_dir, "and figures in", fig_dir, "\n")
