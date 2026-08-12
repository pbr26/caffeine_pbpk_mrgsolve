# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# _setup.R -- shared configuration and helpers. Sourced by every numbered script.
# Assumes the working directory is the project root.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

set.seed(20260806)

PROJ    <- normalizePath(".", mustWork = TRUE)
DIR_MOD <- file.path(PROJ, "models")
DIR_DAT <- file.path(PROJ, "data")
DIR_PHYS<- file.path(DIR_DAT, "physiology")
DIR_DIG <- file.path(DIR_DAT, "digitized")
DIR_FIG <- file.path(PROJ, "figures")
DIR_TAB <- file.path(PROJ, "tables")

for (d in c(DIR_MOD, DIR_PHYS, DIR_DIG, DIR_FIG, DIR_TAB)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

theme_set(theme_bw(base_size = 11))

# --- Helpers -----------------------------------------------------------------

save_fig <- function(plot, name, width = 8, height = 5, dpi = 200) {
  path <- file.path(DIR_FIG, paste0(name, ".png"))
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi)
  message("[fig] ", path)
  invisible(path)
}

save_tab <- function(df, name) {
  path <- file.path(DIR_TAB, paste0(name, ".csv"))
  readr::write_csv(df, path)
  message("[tab] ", path)
  invisible(path)
}

#' Fold-error between prediction and observation (a standard PBPK metric).
#' A value within [0.5, 2] (2-fold) is the usual PBPK acceptance band.
fold_error <- function(pred, obs) pred / obs

#' Absolute average fold error across paired points.
aafe <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs) & obs > 0 & pred > 0
  10 ^ mean(abs(log10(pred[ok] / obs[ok])))
}
