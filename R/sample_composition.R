suppressPackageStartupMessages({ library(tidyverse); library(arrow) })

out_dir <- "results"
theta <- read_parquet(file.path(out_dir, "features_all_pipelines.parquet")) %>%
  filter(band == "theta", outcome %in% c("mean", "sd"))
feat <- theta %>% filter(outcome == "mean")

pipe <- c("reference", "artifact", "epoch_sec", "normalise", "aperiodic", "roi")

diffs <- feat %>%
  select(all_of(pipe), subject, condition, value) %>%
  pivot_wider(names_from = condition, values_from = value) %>%
  filter(is.finite(SD), is.finite(NS)) %>%
  mutate(diff = SD - NS)

curve <- function(d) d %>% group_by(across(all_of(pipe))) %>%
  summarise(n = n(), dz = mean(diff) / sd(diff), .groups = "drop")

share <- function(cv, node) {
  tot <- sum((cv$dz - mean(cv$dz))^2)
  cv %>% group_by(.data[[node]]) %>% summarise(m = mean(dz), k = n(), .groups = "drop") %>%
    summarise(s = 100 * sum(k * (m - mean(cv$dz))^2) / tot) %>% pull(s)
}
nodes <- c("reference", "artifact", "roi", "normalise", "aperiodic", "epoch_sec")

full <- curve(diffs)
full_shares <- sapply(nodes, function(k) share(full, k))
cat("Unrestricted theta-mean multiverse, main-effect shares (%):\n"); print(round(full_shares, 1))

loo <- map_dfr(sort(unique(diffs$subject)), function(s) {
  cv <- curve(filter(diffs, subject != s))
  tibble(dropped = s, n_specs = nrow(cv), median_dz = median(cv$dz),
         r_n_dz = cor(cv$n, cv$dz),
         ica_median = median(cv$dz[cv$artifact == "ica"]),
         threshold_median = median(cv$dz[cv$artifact == "threshold"])) %>%
    bind_cols(as_tibble_row(sapply(nodes, function(k) share(cv, k))))
})
write_csv(loo, file.path(out_dir, "leave_one_out_influence.csv"))
cat("\nLeave-one-participant-out, most influential by reference share:\n")
print(loo %>% arrange(reference) %>% select(dropped, reference, artifact, roi, normalise, r_n_dz, threshold_median, ica_median) %>% head(5), digits = 3)

two_s <- diffs %>% filter(epoch_sec == 2)
present <- two_s %>% distinct(reference, artifact, subject) %>% count(subject) %>% filter(n == 6) %>% pull(subject)
six <- two_s %>% distinct(reference, artifact) %>% arrange(reference, artifact) %>%
  pmap_dfr(function(reference, artifact) {
    d <- two_s[two_s$reference == reference & two_s$artifact == artifact, ]
    own <- curve(d); com <- curve(filter(d, subject %in% present))
    tibble(reference = reference, artifact = artifact,
           n_own = n_distinct(d$subject), n_common = length(present),
           own_median_dz = median(own$dz), common_median_dz = median(com$dz))
  })
write_csv(six, file.path(out_dir, "six_pipelines_own_vs_common.csv"))
cat("\nSix 2-s pipelines, own sample vs the", length(present), "common participants:\n"); print(six, digits = 3)

worst <- loo %>% arrange(reference) %>% slice(1) %>% pull(dropped)
canon <- theta %>% filter(reference == "average", epoch_sec == 2, roi == "centro_temporal",
                          normalise == "absolute", aperiodic == "retained") %>%
  select(subject, artifact, condition, outcome, value, reject_pct) %>%
  pivot_wider(names_from = outcome, values_from = value) %>%
  rename(theta_power = mean, theta_sd = sd)
others <- canon %>% filter(subject != worst) %>% group_by(artifact, condition) %>%
  summarise(n_others = n(), others_median = median(theta_power), others_mad = mad(theta_power),
            others_p5 = quantile(theta_power, 0.05), others_p95 = quantile(theta_power, 0.95), .groups = "drop")
rec <- canon %>% filter(subject == worst) %>%
  left_join(others, by = c("artifact", "condition")) %>%
  mutate(modified_z = (theta_power - others_median) / others_mad) %>%
  select(subject, artifact, condition, theta_power, theta_sd, reject_pct, n_others, others_median, others_mad,
         others_p5, others_p95, modified_z) %>%
  arrange(artifact, condition)
write_csv(rec, file.path(out_dir, "influential_recording.csv"))
cat("\nMost influential participant (", worst, "), centro-temporal absolute theta power, average reference, 2 s.\n",
    "theta_sd is the SD across epochs; the others_ columns describe the remaining participants under the same pipeline.\n", sep = "")
print(rec %>% select(artifact, condition, theta_power, theta_sd, others_median, others_p95, modified_z), digits = 3)
