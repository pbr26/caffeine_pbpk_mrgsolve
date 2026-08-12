# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# The mrgsolve model: compiles, conserves mass when elimination is off, and
# reproduces the validated fit within 2-fold at the calibrated parameters.

test_that("the model compiles and has the expected compartments", {
  m <- .build_model()
  cmts <- mrgsolve::cmt(m)
  expect_true(all(c("DEPOT","LI","GU","AR","VE","MET","URINE") %in% cmts))
})

test_that("mass is conserved when metabolism and renal CL are switched off", {
  m <- .build_model() %>% mrgsolve::param(Vmax = 0, CLr = 0)
  out <- m %>% mrgsolve::ev(amt = 250, cmt = "DEPOT") %>%
    mrgsolve::mrgsim(end = 24, delta = 1) %>% as.data.frame()
  amt_cols <- c("DEPOT","AD","BO","BR","HE","KI","LI","LU","MU","SK","TH",
                "GU","RE","AR","VE")
  total <- rowSums(out[, amt_cols])
  expect_true(all(abs(total - 250) < 250 * 1e-3))   # <0.1% drift over 24 h
})

test_that("with elimination on, drug is cleared into the metabolism sink", {
  m <- .build_model()
  out <- m %>% mrgsolve::ev(amt = 250, cmt = "DEPOT") %>%
    mrgsolve::mrgsim(end = 48, delta = 1) %>% as.data.frame()
  expect_gt(tail(out$MET, 1), 0)                    # some drug metabolized
  expect_lt(tail(out$CP, 1), max(out$CP))           # concentration declines
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
  pred <- approx(s$time, s$CP, xout = obs$time_h)$y
  aafe <- 10 ^ mean(abs(log10(pred / obs$conc_ug_mL)), na.rm = TRUE)
  expect_lt(aafe, 2)
})
