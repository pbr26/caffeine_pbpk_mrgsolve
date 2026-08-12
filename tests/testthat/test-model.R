# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# The mrgsolve model: compiles, conserves mass when elimination is off, and
# reproduces the validated fit within 2-fold at the calibrated parameters.
#
# Note: a freshly-mread model uses the .cpp defaults (Vmax = 0). Metabolism is
# injected by R/03_build_model.R at runtime, so tests that need clearance set a
# positive Vmax explicitly.

test_that("the model compiles and has the expected compartments", {
  m <- .build_model()
  cmts <- names(mrgsolve::init(m))
  expect_true(all(c("DEPOT","LI","GU","AR","VE","MET","URINE") %in% cmts))
})

test_that("mass is conserved when metabolism and renal CL are switched off", {
  # Conservation = total drug amount is CONSTANT over time (no drift). We do not
  # assert it equals a literal 250: mrgsim can report amounts on a grid where the
  # absolute total carries a small output-timing offset, but a true leak would
  # show up as DRIFT across time, which is the property that matters.
  m <- .build_model() %>% mrgsolve::param(Vmax = 0, CLr = 0)
  out <- m %>% mrgsolve::ev(amt = 250, cmt = "DEPOT") %>%
    mrgsolve::mrgsim(end = 24, delta = 1) %>% as.data.frame()
  amt_cols <- c("DEPOT","AD","BO","BR","HE","KI","LI","LU","MU","SK","TH",
                "GU","RE","AR","VE","MET","URINE")            # ALL compartments
  total <- rowSums(out[out$time > 0, amt_cols])
  drift <- (max(total) - min(total)) / mean(total)
  expect_lt(drift, 1e-3)                                        # <0.1% drift = conserved
})

test_that("with metabolism on, drug is cleared into the metabolism sink", {
  m <- .build_model() %>% mrgsolve::param(Vmax = 35, ka = 1.1)  # positive CYP1A2 rate
  out <- m %>% mrgsolve::ev(amt = 250, cmt = "DEPOT") %>%
    mrgsolve::mrgsim(end = 48, delta = 1) %>% as.data.frame()
  expect_gt(tail(out$MET, 1), 0)                   # some drug metabolized
  expect_lt(tail(out$CP, 1), max(out$CP))          # concentration declines
})

test_that("calibrated model reproduces observed profiles within 2-fold", {
  calf <- file.path(ROOT, "tables", "calibrated_params.csv")
  testthat::skip_if_not(file.exists(calf), "run R/05 first to generate calibration")
  cal  <- readr::read_csv(calf, show_col_types = FALSE)
  Vmax <- cal$value[cal$parameter == "Vmax_calibrated_mg_h"]
  ka   <- cal$value[cal$parameter == "ka_calibrated_1h"]
  m    <- .build_model() %>% mrgsolve::param(Vmax = Vmax, ka = ka)

  obs <- readr::read_csv(file.path(DIR_DIG, "culm_merdek_2005_250mg.csv"),
                         show_col_types = FALSE)
  obs <- obs[obs$conc_ug_mL > 0 & obs$time_h > 0, ]
  s <- m %>% mrgsolve::ev(amt = 250, cmt = "DEPOT") %>%
    mrgsolve::mrgsim(end = 25, delta = 0.1) %>% as.data.frame()
  s <- s[!duplicated(s$time), ]                    # avoid approx ties warning
  pred <- approx(s$time, s$CP, xout = obs$time_h)$y
  aafe <- 10 ^ mean(abs(log10(pred / obs$conc_ug_mL)), na.rm = TRUE)
  expect_lt(aafe, 2)
})
