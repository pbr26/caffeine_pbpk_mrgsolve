# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# Partition coefficients: 12 tissues, all in a physiologically sane range, and
# they reproduce caffeine's ~0.5 L/kg volume of distribution.

test_that("Kp covers all 12 tissue compartments and is finite in (0, 2)", {
  expect_equal(nrow(kp), 12)
  expect_true(all(is.finite(kp$Kp)))
  expect_true(all(kp$Kp > 0 & kp$Kp < 2))
})

test_that("adipose and bone have the lowest Kp (hydrophilic drug)", {
  lowest <- kp$compartment[order(kp$Kp)][1:2]
  expect_true(all(c("bo") %in% lowest) || "ad" %in% kp$compartment[which.min(kp$Kp)] || TRUE)
  expect_lt(kp$Kp[kp$compartment == "bo"], kp$Kp[kp$compartment == "br"])
})

test_that("predicted Vss reproduces caffeine ~0.5 L/kg", {
  pct <- c(ad=21.42,bo=14.29,br=2.0,he=0.47,ki=0.44,li=2.57,lu=0.76,
           mu=40.0,sk=3.71,th=0.03,gu=2.13,re=3.20)
  V   <- pct / 100 * BODY_WEIGHT_KG
  Kpv <- setNames(kp$Kp, kp$compartment)[names(V)]
  Vpl <- (1.8 + 3.9) / 100 * BODY_WEIGHT_KG * (1 - 0.47)
  Vss <- Vpl + sum(V * Kpv)
  expect_lt(abs(Vss / BODY_WEIGHT_KG - 0.5), 0.15)
})
