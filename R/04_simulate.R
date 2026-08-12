# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# 04_simulate.R -- simulate single oral caffeine doses and overlay the
# venous-plasma prediction on the observed profiles (digitized from OSP).
#
# Depends on 03_build_model.R having compiled `mod`. Run in the same session,
# or this script rebuilds it via source().
#
# Doses simulated (to match the two observed studies):
#   250 mg single oral -> Culm-Merdek et al. 2005
#   200 mg single oral -> Jeppesen et al. 1996 (EM)
# Output: figures/sim_vs_observed.png and tables/sim_predictions.csv
# -----------------------------------------------------------------------------

source("R/_setup.R")
if (!exists("mod")) source("R/03_build_model.R")

# --- Observed profiles -------------------------------------------------------
obs_files <- c(
  "250" = file.path(DIR_DIG, "culm_merdek_2005_250mg.csv"),
  "200" = file.path(DIR_DIG, "jeppesen_1996_200mg.csv")
)
obs <- dplyr::bind_rows(lapply(names(obs_files), function(d) {
  readr::read_csv(obs_files[[d]], show_col_types = FALSE) %>%
    dplyr::transmute(dose_mg = as.numeric(d), time_h, conc = conc_ug_mL,
                     sd = sd_ng_mL / 1000, study)   # sd ng/mL -> ug/mL
}))

# --- Simulate each dose ------------------------------------------------------
sim_one <- function(dose_mg, tend = 25) {
  d <- ev(amt = dose_mg, cmt = "DEPOT")
  mod %>%
    ev(d) %>%
    mrgsim(end = tend, delta = 0.05) %>%
    as.data.frame() %>%
    dplyr::transmute(dose_mg = dose_mg, time_h = time, pred = CP)
}
sim <- dplyr::bind_rows(lapply(c(250, 200), sim_one))
save_tab(sim, "sim_predictions")

# --- Overlay figure ----------------------------------------------------------
lab <- c("250" = "250 mg (Culm-Merdek 2005)", "200" = "200 mg (Jeppesen 1996)")
sim$grp <- lab[as.character(sim$dose_mg)]
obs$grp <- lab[as.character(obs$dose_mg)]

p <- ggplot() +
  geom_line(data = sim, aes(time_h, pred, colour = grp), linewidth = 0.8) +
  geom_point(data = obs, aes(time_h, conc, colour = grp), size = 2) +
  geom_errorbar(data = obs, aes(time_h, ymin = pmax(conc - sd, 0),
                                ymax = conc + sd, colour = grp),
                width = 0.3, na.rm = TRUE) +
  labs(x = "Time after dose (h)", y = "Caffeine plasma conc (ug/mL)",
       colour = NULL,
       title = "Caffeine PBPK: venous-plasma prediction vs observed",
       subtitle = "Perfusion-limited whole-body model (mrgsolve); Kp/ka provisional") +
  theme(legend.position = "top")
save_fig(p, "sim_vs_observed", width = 8, height = 5)

# --- Quick fit preview (formal validation is step 5) -------------------------
preview <- obs %>%
  dplyr::group_by(dose_mg) %>%
  dplyr::group_modify(function(df, key) {
    s <- sim[sim$dose_mg == key$dose_mg, ]
    pr <- approx(s$time_h, s$pred, xout = df$time_h)$y
    tibble::tibble(n = sum(is.finite(pr) & df$conc > 0),
                   AAFE = aafe(pr, df$conc),
                   Cmax_obs = max(df$conc),
                   Cmax_pred = max(s$pred))
  }) %>% dplyr::ungroup()
save_tab(preview, "sim_fit_preview")
print(preview)

message("04_simulate.R complete. See figures/sim_vs_observed.png")
