# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# Caffeine compound parameters: physical plausibility + provenance discipline.

gv <- function(p) suppressWarnings(as.numeric(cp$value[cp$parameter == p]))

test_that("molecular weight matches caffeine (194.19 g/mol)", {
  expect_equal(gv("molecular_weight"), 194.19, tolerance = 1)
})

test_that("fraction unbound is a valid fraction", {
  fu <- gv("fraction_unbound_plasma")
  expect_true(fu > 0 && fu <= 1)
})

test_that("CYP1A2 kinetics are positive", {
  expect_gt(gv("CYP1A2_Km"), 0)
  expect_gt(gv("CYP1A2_Vmax_invitro"), 0)
})

test_that("every parameter carries a status and a source", {
  expect_true(all(nzchar(cp$status)))
  expect_true(all(nzchar(cp$source)))
  expect_true(all(cp$status %in%
    c("verified","confirm-units","confirm-alternative","todo","derived")))
})
