# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# Physiology provenance + balance: a table that fails these is a data bug.

test_that("physiology has the expected 14 compartments", {
  expect_equal(nrow(phys), 14)
  expect_true(all(c("adipose","muscle","liver","lung","gut",
                    "arterial_blood","venous_blood") %in% phys$compartment))
})

test_that("compartment volumes sum to ~ body weight", {
  V <- phys$pct_body_weight / 100 * BODY_WEIGHT_KG
  expect_lt(abs(sum(V) - BODY_WEIGHT_KG) / BODY_WEIGHT_KG, 0.10)
})

test_that("arterial inflows sum to ~ cardiac output (lung + blood excluded)", {
  art <- phys[!phys$compartment %in% c("lung","arterial_blood","venous_blood"), ]
  Q <- art$pct_cardiac_output / 100 * CARDIAC_OUTPUT_LH
  expect_lt(abs(sum(Q) - CARDIAC_OUTPUT_LH) / CARDIAC_OUTPUT_LH, 0.10)
})

test_that("liver flow is hepatic-artery-only and gut carries the portal share", {
  # documented convention: liver row = HA (4.6%), gut row = portal (18.1%)
  expect_equal(phys$pct_cardiac_output[phys$compartment == "liver"], 4.6, tolerance = 0.1)
  expect_equal(phys$pct_cardiac_output[phys$compartment == "gut"], 18.1, tolerance = 0.1)
})
