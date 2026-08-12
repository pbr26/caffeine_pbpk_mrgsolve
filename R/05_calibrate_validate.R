# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# 05_calibrate_validate.R -- calibrate the CYP1A2 catalytic rate (and oral ka)
# to the observed caffeine profiles, then validate.
#
# Rationale / faithfulness: distribution is already fixed mechanistically (PK-Sim
# Schmitt Kp -> Vss ~0.5 L/kg). The remaining misfit is a pure terminal-slope
# (clearance) problem: the in-vitro CYP1A2 Vmax (flagged 'confirm-units') gives a
# ~1.8 h half-life vs caffeine's ~5 h. Britz et al. 2019 themselves state they
# "fit the catalytic rate constants of ... CYP1A2" to clinical data -- so we do
# the same: optimise a scale factor on Vmax plus the absorption rate ka, jointly
# across both observed studies, in log space.
#
# Fixed (NOT fitted): all physiology, Kp, fu, Km, MW -- provenance preserved.
# Fitted: kcat_scale (multiplies whole-liver Vmax), ka. Two parameters, ~23 obs.
#
# Depends on a compiled `mod` from 03_build_model.R.
# Outputs: tables/calibrated_params.csv, tables/validation_metrics.csv,
#          figures/calibrated_fit.png, figures/gof.png
# -----------------------------------------------------------------------------

source("R/_setup.R")
if (!exists("mod")) source("R/03_build_model.R")

# --- Observed data -----------------------------------------------------------
obs <- dplyr::bind_rows(
  readr::read_csv(file.path(DIR_DIG, "culm_merdek_2005_250mg.csv"), show_col_types = FALSE) %>%
    dplyr::transmute(dose = 250, time_h, conc = conc_ug_mL),
  readr::read_csv(file.path(DIR_DIG, "jeppesen_1996_200mg.csv"), show_col_types = FALSE) %>%
    dplyr::transmute(dose = 200, time_h, conc = conc_ug_mL)
) %>% dplyr::filter(conc > 0, time_h > 0)        # drop t=0 / non-positive

Vmax0 <- as.numeric(param(mod)$Vmax)             # in-vitro-scaled starting Vmax
ka0   <- as.numeric(param(mod)$ka)

# --- Simulate one dose at a given (kcat_scale, ka) ---------------------------
sim_dose <- function(dose, kscale, ka, times) {
  mod %>%
    param(Vmax = Vmax0 * kscale, ka = ka) %>%
    ev(amt = dose, cmt = "DEPOT") %>%
    mrgsim(end = max(times) + 1, delta = 0.1) %>%
    as.data.frame() -> s
  approx(s$time, s$CP, xout = times)$y
}

# --- Objective: joint weighted least squares in log space --------------------
obj <- function(logpar) {
  kscale <- exp(logpar[1]); ka <- exp(logpar[2])
  ss <- 0
  for (d in unique(obs$dose)) {
    o <- obs[obs$dose == d, ]
    p <- sim_dose(d, kscale, ka, o$time_h)
    ok <- is.finite(p) & p > 0
    ss <- ss + sum((log(p[ok]) - log(o$conc[ok]))^2)
  }
  ss
}

# --- Optimise (start: reduce Vmax ~3x to hit ~5 h half-life) -----------------
fit <- optim(c(log(0.35), log(ka0)), obj, method = "Nelder-Mead",
             control = list(reltol = 1e-8, maxit = 500))
kscale_hat <- exp(fit$par[1]); ka_hat <- exp(fit$par[2])
Vmax_hat   <- Vmax0 * kscale_hat

# --- Fitted disposition summary ----------------------------------------------
modC <- mod %>% param(Vmax = Vmax_hat, ka = ka_hat)
# terminal half-life from a fine single-dose sim
fine <- modC %>% ev(amt = 250, cmt = "DEPOT") %>% mrgsim(end = 36, delta = 0.1) %>% as.data.frame()
tail <- fine[fine$time >= 12 & fine$time <= 30, ]
lam  <- -coef(lm(log(CP) ~ time, tail))[2]
thalf <- as.numeric(log(2) / lam)

cal <- tibble::tibble(
  parameter = c("kcat_scale_on_Vmax", "Vmax_calibrated_mg_h", "ka_calibrated_1h",
                "terminal_half_life_h"),
  value  = c(kscale_hat, Vmax_hat, ka_hat, thalf),
  note   = c("multiplier on in-vitro-scaled whole-liver Vmax (fitted, as Britz fit kcat)",
             sprintf("was %.1f (in-vitro scaled); calibrated to clinical data", Vmax0),
             sprintf("was %.2f (provisional)", ka0),
             "caffeine literature ~4-6 h"))
save_tab(cal, "calibrated_params")

# --- Validation metrics (fold error, AAFE, % within 2-fold) ------------------
val <- obs %>% dplyr::group_by(dose) %>% dplyr::group_modify(function(df, key) {
  p  <- sim_dose(key$dose, kscale_hat, ka_hat, df$time_h)
  fe <- p / df$conc
  tibble::tibble(n = sum(is.finite(fe)),
                 AAFE = aafe(p, df$conc),
                 pct_within_2fold = 100 * mean(fe >= 0.5 & fe <= 2, na.rm = TRUE),
                 Cmax_obs = max(df$conc), Cmax_pred = max(p, na.rm = TRUE))
}) %>% dplyr::ungroup()
save_tab(val, "validation_metrics")
print(cal); print(val)

# --- Figures: calibrated overlay + goodness-of-fit ---------------------------
lab <- c("250" = "250 mg (Culm-Merdek 2005)", "200" = "200 mg (Jeppesen 1996)")
curve <- dplyr::bind_rows(lapply(unique(obs$dose), function(d) {
  s <- modC %>% ev(amt = d, cmt = "DEPOT") %>% mrgsim(end = 25, delta = 0.1) %>% as.data.frame()
  tibble::tibble(grp = lab[as.character(d)], time_h = s$time, pred = s$CP)
}))
obs$grp <- lab[as.character(obs$dose)]
p1 <- ggplot() +
  geom_line(data = curve, aes(time_h, pred, colour = grp), linewidth = 0.8) +
  geom_point(data = obs, aes(time_h, conc, colour = grp), size = 2) +
  labs(x = "Time after dose (h)", y = "Caffeine plasma (ug/mL)", colour = NULL,
       title = "Caffeine PBPK - calibrated fit vs observed",
       subtitle = sprintf("Schmitt Kp (Vss~0.5 L/kg); CYP1A2 rate fitted; t1/2=%.1f h", thalf)) +
  theme(legend.position = "top")
save_fig(p1, "calibrated_fit")

gof <- obs %>% dplyr::group_by(dose) %>% dplyr::group_modify(function(df, key) {
  df$pred <- sim_dose(key$dose, kscale_hat, ka_hat, df$time_h); df
}) %>% dplyr::ungroup()
rng <- range(c(gof$conc, gof$pred), na.rm = TRUE)
p2 <- ggplot(gof, aes(conc, pred, colour = lab[as.character(dose)])) +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = 1, intercept = c(log10(2), -log10(2)), linetype = 2) +
  geom_point(size = 2) + scale_x_log10() + scale_y_log10() +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = "Observed (ug/mL)", y = "Predicted (ug/mL)", colour = NULL,
       title = "Goodness of fit (dashed = 2-fold)") + theme(legend.position = "top")
save_fig(p2, "gof")

message(sprintf("05 complete. Calibrated: Vmax %.1f->%.1f mg/h, ka %.2f->%.2f /h, t1/2 %.1f h.",
                Vmax0, Vmax_hat, ka0, ka_hat, thalf))
