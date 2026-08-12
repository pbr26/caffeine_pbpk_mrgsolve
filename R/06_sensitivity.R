# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# 06_sensitivity.R -- local (one-at-a-time) sensitivity analysis of the
# calibrated caffeine PBPK model.
#
# For each key parameter we perturb it by +/- DELTA and record the normalized
# sensitivity coefficient  S = (dM/M) / (dP/P)  for three PK metrics of a 250 mg
# oral dose: Cmax, AUC(0-48h), and terminal half-life. |S| ~ 1 means a metric
# scales roughly proportionally with the parameter; |S| ~ 0 means insensitive.
#
# Baseline = the model with the step-05 calibrated Vmax and ka.
# Outputs: tables/sensitivity.csv, figures/tornado_cmax.png, figures/tornado_auc.png
# -----------------------------------------------------------------------------

source("R/_setup.R")
if (!exists("mod")) source("R/03_build_model.R")

# --- Apply the calibrated baseline (from step 05) ----------------------------
calf <- file.path(DIR_TAB, "calibrated_params.csv")
if (file.exists(calf)) {
  cal <- readr::read_csv(calf, show_col_types = FALSE)
  Vmax_base <- cal$value[cal$parameter == "Vmax_calibrated_mg_h"]
  ka_base   <- cal$value[cal$parameter == "ka_calibrated_1h"]
  mod <- mod %>% param(Vmax = Vmax_base, ka = ka_base)
  message("Sensitivity around calibrated baseline (Vmax=", round(Vmax_base,1),
          ", ka=", round(ka_base,2), ").")
} else {
  message("No calibrated_params.csv -> using build defaults as baseline.")
}

DOSE <- 250; DELTA <- 0.15

# --- PK metrics for a single 250 mg oral dose --------------------------------
metrics <- function(m) {
  s <- m %>% ev(amt = DOSE, cmt = "DEPOT") %>% mrgsim(end = 48, delta = 0.05) %>% as.data.frame()
  auc  <- sum(diff(s$time) * (head(s$CP,-1) + tail(s$CP,-1)) / 2)   # trapezoid 0-48h
  tl   <- s[s$time >= 12 & s$time <= 36 & s$CP > 0, ]              # positive tail only
  lam  <- -coef(lm(log(CP) ~ time, tl))[2]
  c(Cmax = max(s$CP), AUC = auc, thalf = as.numeric(log(2)/lam))
}

# robust named-vector view of the current parameters
P0     <- as.list(param(mod))
# Only parameters that keep the model valid one-at-a-time. Individual branch
# flows are NOT perturbed alone (they must sum to cardiac output, or blood-flow
# conservation breaks); instead Q_ALL scales the whole circulation together.
scalar <- c("ka","Vmax","Km","fup","BP")                        # multiply by (1+/-DELTA)
kp_all <- paste0("Kp_", c("ad","bo","br","he","ki","li","lu","mu","sk","th","gu","re"))
q_all  <- c("Q_ad","Q_bo","Q_br","Q_he","Q_ki","Q_ha","Q_mu","Q_sk","Q_th","Q_gu","Q_re","Q_co")
base   <- metrics(mod)

perturb <- function(par, factor) {
  grp <- switch(par, "Kp_ALL" = kp_all, "Q_ALL" = q_all, NULL)
  if (!is.null(grp)) {
    new <- setNames(vapply(grp, function(k) as.numeric(P0[[k]]) * factor, numeric(1)), grp)
  } else {
    new <- setNames(as.numeric(P0[[par]]) * factor, par)
  }
  param(mod, new)
}

probes <- c(scalar, "Kp_ALL", "Q_ALL")
rows <- lapply(probes, function(par) {
  hi <- metrics(perturb(par, 1 + DELTA))
  lo <- metrics(perturb(par, 1 - DELTA))
  # central normalized sensitivity for each metric
  s <- (hi - lo) / base / (2 * DELTA)
  tibble::tibble(parameter = par, metric = names(base),
                 S = as.numeric(s), base = as.numeric(base),
                 hi = as.numeric(hi), lo = as.numeric(lo))
})
sens <- dplyr::bind_rows(rows)
save_tab(sens, "sensitivity")

# --- Tornado plots for Cmax and AUC ------------------------------------------
tornado <- function(metric_name, fname, ttl) {
  d <- sens %>% dplyr::filter(metric == metric_name) %>%
    dplyr::mutate(parameter = reorder(parameter, abs(S)))
  p <- ggplot(d, aes(parameter, S, fill = S > 0)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#2c7fb8", "FALSE" = "#de2d26"), guide = "none") +
    labs(x = NULL, y = "Normalized sensitivity  S = (dM/M)/(dP/P)",
         title = ttl, subtitle = "250 mg oral; +/-15% one-at-a-time") +
    geom_hline(yintercept = 0, linewidth = 0.3)
  save_fig(p, fname, width = 7, height = 4.5)
}
tornado("Cmax", "tornado_cmax", "Sensitivity of caffeine Cmax")
tornado("AUC",  "tornado_auc",  "Sensitivity of caffeine AUC(0-48h)")

print(sens %>% dplyr::filter(metric %in% c("Cmax","AUC")) %>%
        dplyr::arrange(metric, -abs(S)))
message("06_sensitivity.R complete. Tables/tornado figures written.")
