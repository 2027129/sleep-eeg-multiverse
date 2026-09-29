data_dir      <- "ds004902"
N_CORES       <- 4

FILTER_HP     <- 0.5
FILTER_LP     <- 45
ARTIFACT_UV   <- 100
BAD_CHAN_Z    <- 3

FORK_REFERENCE  <- c("average", "mastoid", "cz")
FORK_ARTIFACT   <- c("threshold", "ica")
FORK_EPOCH_SEC  <- c(2, 4, 8, 10)
FORK_NORMALISE  <- c("absolute", "rel_1_30", "rel_4_30")
FORK_APERIODIC  <- c("retained", "removed")
FORK_OUTCOME    <- c("mean", "sd", "cv")

BANDS <- list(theta = c(4, 8), alpha = c(8, 13), beta = c(13, 30))

ROI_CUI <- list(
  frontal = c("Fp1","Fpz","Fp2","AF7","AF3","AF4","AF8",
              "F7","F5","F3","F1","Fz","F2","F4","F6","F8"),
  centro_temporal = c("FT7","FC5","FC3","FC1","FC2","FC4","FC6","FT8",
                      "T7","C5","C3","C1","Cz","C2","C4","C6","T8",
                      "TP7","CP5","CP3","CP1","CPz","CP2","CP4","CP6",
                      "TP8","TP9","TP10"),
  parieto_occipital = c("P7","P5","P3","P1","Pz","P2","P4","P6","P8",
                        "PO7","PO3","POz","PO4","PO8","O1","Oz","O2")
)

pkgs <- c("tidyverse","fs","jsonlite","signal","fastICA",
          "furrr","future","arrow")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(tidyverse); library(fs); library(jsonlite); library(furrr)

output_dir <- "results"
cache_dir  <- file.path(output_dir, "cache")
dir_create(output_dir); dir_create(cache_dir)

read_bids_eeg <- function(set_path) {
  base <- sub("_eeg\\.set$", "", set_path)
  jf <- paste0(base, "_eeg.json")
  cf <- paste0(base, "_channels.tsv")
  ff <- paste0(base, "_eeg.fdt")
  if (!all(file.exists(jf, cf, ff))) return(NULL)

  meta  <- fromJSON(jf)
  chans <- suppressMessages(read_tsv(cf, show_col_types = FALSE))
  fs_hz  <- meta$SamplingFrequency
  nbchan <- nrow(chans)

  pnts <- (file.size(ff) %/% 4) %/% nbchan
  if (pnts < fs_hz * 30) return(NULL)

  raw <- readBin(ff, "numeric", size = 4, n = pnts * nbchan, endian = "little")
  mat <- matrix(raw, nrow = nbchan, ncol = pnts)
  rownames(mat) <- chans$name

  keep <- if ("type" %in% names(chans)) which(toupper(chans$type) == "EEG") else seq_len(nbchan)
  list(data = mat[keep, , drop = FALSE], fs = fs_hz,
       names = chans$name[keep], pnts = pnts)
}

bandpass <- function(x, fs) {
  hp <- signal::butter(2, FILTER_HP / (fs/2), type = "high")
  lp <- signal::butter(4, FILTER_LP / (fs/2), type = "low")
  tryCatch(signal::filtfilt(lp, signal::filtfilt(hp, x)),
           error = function(e) rep(NA_real_, length(x)))
}

detect_bad <- function(mat) {
  lv <- log(apply(mat, 1, var, na.rm = TRUE) + 1e-12)
  m <- median(lv); s <- mad(lv)
  if (!is.finite(s) || s == 0) return(character(0))
  rownames(mat)[abs((lv - m) / s) > BAD_CHAN_Z]
}

apply_reference <- function(mat, scheme) {
  if (scheme == "average") {
    sweep(mat, 2, colMeans(mat), "-")
  } else if (scheme == "mastoid") {
    mast <- intersect(c("TP9","TP10"), rownames(mat))
    if (length(mast) == 0) return(NULL)
    sweep(mat, 2, colMeans(mat[mast, , drop = FALSE]), "-")
  } else if (scheme == "cz") {
    if (!"Cz" %in% rownames(mat)) return(NULL)
    sweep(mat, 2, mat["Cz", ], "-")
  }
}

remove_ocular_ica <- function(mat, n_comp = 20, thresh = 0.5) {
  eog_ch <- intersect(c("Fp1","Fp2"), rownames(mat))
  if (length(eog_ch) == 0) return(mat)
  eog <- colMeans(mat[eog_ch, , drop = FALSE])

  fit <- tryCatch(
    fastICA::fastICA(t(mat), n.comp = min(n_comp, nrow(mat)), method = "C"),
    error = function(e) NULL)
  if (is.null(fit)) return(mat)

  r <- apply(fit$S, 2, function(s) abs(cor(s, eog)))
  drop <- which(r > thresh)
  if (length(drop) == 0) return(mat)

  S <- fit$S; S[, drop] <- 0
  cleaned <- t(S %*% fit$A)
  rownames(cleaned) <- rownames(mat)
  cleaned
}

stage_a <- function(set_path) {
  eeg <- read_bids_eeg(set_path)
  if (is.null(eeg)) return(NULL)

  filt <- t(apply(eeg$data, 1, bandpass, fs = eeg$fs))
  rownames(filt) <- eeg$names
  filt <- filt[complete.cases(filt), , drop = FALSE]
  if (nrow(filt) < 20) return(NULL)

  good <- setdiff(rownames(filt), detect_bad(filt))
  filt <- filt[good, , drop = FALSE]

  out <- list()
  for (ref in FORK_REFERENCE) {
    reref <- apply_reference(filt, ref)
    if (is.null(reref)) next
    for (art in FORK_ARTIFACT) {
      cleaned <- if (art == "ica") remove_ocular_ica(reref) else reref
      out[[paste(ref, art, sep = "|")]] <- cleaned
    }
  }
  list(variants = out, fs = eeg$fs)
}

epoch_psd <- function(mat, fs, epoch_sec) {
  L <- round(epoch_sec * fs)
  starts <- seq(1, ncol(mat) - L + 1, by = L)
  win <- 0.5 - 0.5 * cos(2*pi*(0:(L-1))/(L-1))
  wn  <- sum(win^2)

  freqs <- (0:(L-1)) * fs / L
  keep  <- which(freqs >= 1 & freqs <= 45)

  roi_names <- names(ROI_CUI)
  roi_idx <- lapply(ROI_CUI, function(chs) which(rownames(mat) %in% chs))

  ok <- c()
  acc <- array(0, dim = c(length(starts), length(roi_names), length(keep)))

  for (i in seq_along(starts)) {
    seg <- mat[, starts[i]:(starts[i]+L-1), drop = FALSE]
    if (max(abs(seg)) > ARTIFACT_UV) next
    ok <- c(ok, i)
    for (r in seq_along(roi_names)) {
      idx <- roi_idx[[r]]
      if (length(idx) == 0) next
      p <- sapply(idx, function(ch) {
        (Mod(fft(seg[ch, ] * win))^2)[keep] / (fs * wn)
      })
      acc[i, r, ] <- rowMeans(p) * 2
    }
  }

  if (length(ok) < 10) return(NULL)
  list(psd = acc[ok, , , drop = FALSE], freqs = freqs[keep],
       n_epochs = length(ok), n_total = length(starts),
       rois = roi_names)
}

band_integral <- function(psd_vec, freqs, lo, hi) {
  idx <- which(freqs >= lo & freqs < hi)
  sum(psd_vec[idx]) * (freqs[2] - freqs[1])
}

remove_aperiodic <- function(psd_vec, freqs, fit_range = c(2, 30)) {
  sel <- which(freqs >= fit_range[1] & freqs <= fit_range[2])
  lf <- log10(freqs[sel]); lp <- log10(psd_vec[sel] + 1e-15)
  cf <- stats::lm.fit(cbind(1, lf), lp)$coefficients
  resid <- rep(0, length(psd_vec))
  resid[sel] <- lp - (cf[1] + cf[2] * lf)
  pmax(resid, 0)
}

stage_c <- function(B) {
  freqs <- B$freqs
  n_ep <- dim(B$psd)[1]
  rows <- list()

  for (aper in FORK_APERIODIC) {
    spec <- B$psd
    if (aper == "removed") {
      for (e in seq_len(n_ep)) for (r in seq_along(B$rois)) {
        spec[e, r, ] <- remove_aperiodic(B$psd[e, r, ], freqs)
      }
    }

    for (norm in FORK_NORMALISE) {
      for (bn in names(BANDS)) {
        lo <- BANDS[[bn]][1]; hi <- BANDS[[bn]][2]

        for (r in seq_along(B$rois)) {
          vals <- vapply(seq_len(n_ep), function(e) {
            num <- band_integral(spec[e, r, ], freqs, lo, hi)
            den <- switch(norm,
              absolute = 1,
              rel_1_30 = band_integral(spec[e, r, ], freqs, 1, 30),
              rel_4_30 = band_integral(spec[e, r, ], freqs, 4, 30))
            if (den == 0) NA_real_ else num / den
          }, numeric(1))

          vals <- vals[is.finite(vals)]
          if (length(vals) < 5) next

          for (outc in FORK_OUTCOME) {
            v <- switch(outc,
              mean = mean(vals),
              sd   = sd(vals),
              cv   = sd(vals) / mean(vals))
            rows[[length(rows)+1]] <- tibble(
              aperiodic = aper, normalise = norm, outcome = outc,
              band = bn, roi = B$rois[r], value = v)
          }
        }
      }
    }
  }
  bind_rows(rows)
}

all_files <- dir_ls(path_expand(data_dir), recurse = TRUE, type = "file")
set_files <- all_files[str_detect(all_files, "task-eyesopen.*_eeg\\.set$")]
if (length(set_files) == 0)
  set_files <- all_files[str_detect(all_files, "_eeg\\.set$")]

cat("Recordings to process:", length(set_files), "\n")
cat("Pipelines per recording:",
    length(FORK_REFERENCE) * length(FORK_ARTIFACT) * length(FORK_EPOCH_SEC) *
    length(FORK_NORMALISE) * length(FORK_APERIODIC) * length(FORK_OUTCOME), "\n\n")

plan(multisession, workers = N_CORES)

process_one <- function(f) {
  cache_f <- file.path(cache_dir, paste0(tools::file_path_sans_ext(basename(f)), ".rds"))
  if (file.exists(cache_f)) return(readRDS(cache_f))

  A <- stage_a(f)
  if (is.null(A)) return(NULL)

  subject <- str_remove(str_extract(basename(f), "sub-[A-Za-z0-9]+"), "sub-")
  session <- str_remove(str_extract(basename(f), "ses-[A-Za-z0-9]+"), "ses-")

  out <- list()
  for (vn in names(A$variants)) {
    parts <- str_split(vn, "\\|")[[1]]
    for (ep in FORK_EPOCH_SEC) {
      B <- epoch_psd(A$variants[[vn]], A$fs, ep)
      if (is.null(B)) next
      feats <- stage_c(B)
      feats$reference <- parts[1]; feats$artifact <- parts[2]
      feats$epoch_sec <- ep
      feats$subject <- subject
      feats$condition <- if (session == "1") "NS" else "SD"
      feats$reject_pct <- 100 * (1 - B$n_epochs / B$n_total)
      out[[length(out)+1]] <- feats
    }
  }
  res <- bind_rows(out)
  saveRDS(res, cache_f)
  res
}

t0 <- Sys.time()
features <- future_map_dfr(set_files, process_one, .progress = TRUE,
                            .options = furrr_options(seed = TRUE))
cat("\nElapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

arrow::write_parquet(features, file.path(output_dir, "features_all_pipelines.parquet"))

spec_id <- c("reference","artifact","epoch_sec","normalise",
             "aperiodic","outcome","band","roi")

paired_effect <- function(d) {
  w <- d %>% select(subject, condition, value) %>%
    pivot_wider(names_from = condition, values_from = value) %>%
    filter(is.finite(NS), is.finite(SD))
  if (nrow(w) < 10) return(tibble(n = nrow(w), dz = NA, p = NA))
  diff <- w$SD - w$NS
  tt <- t.test(diff)
  tibble(n = nrow(w), dz = mean(diff)/sd(diff), p = tt$p.value)
}

spec_curve <- features %>%
  group_by(across(all_of(spec_id))) %>%
  group_modify(~ paired_effect(.x)) %>%
  ungroup() %>%
  mutate(p_fdr = p.adjust(p, "BH"), sig = p_fdr < 0.05)

write_csv(spec_curve, file.path(output_dir, "specification_curve.csv"))

variance_decomp <- spec_curve %>%
  filter(is.finite(dz)) %>%
  { aov(dz ~ reference + artifact + factor(epoch_sec) + normalise +
              aperiodic + outcome + band + roi, data = .) } %>%
  broom::tidy() %>%
  mutate(pct_variance = 100 * sumsq / sum(sumsq)) %>%
  arrange(desc(pct_variance))

write_csv(variance_decomp, file.path(output_dir, "variance_decomposition.csv"))

cat("Specifications run:", nrow(spec_curve), "\n")
cat("Significant after FDR:", sum(spec_curve$sig, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100*mean(spec_curve$sig, na.rm = TRUE)))
cat("\nVariance in effect size explained by each fork:\n")
print(variance_decomp %>% select(term, pct_variance))
cat("\nCui et al.'s own specification:\n")
print(spec_curve %>% filter(reference == "average", artifact == "ica",
                             epoch_sec == 4, normalise == "rel_4_30",
                             aperiodic == "retained", band == "theta"))
