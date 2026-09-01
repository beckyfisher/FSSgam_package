# Golden-master capture for the phase 2 correlation changes.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
out.file <- commandArgs(trailingOnly = TRUE)[1]
set.seed(7)
d <- FSSgam::case_study1
d$site <- factor(d$site)
d$ZONE2 <- factor(ifelse(d$SCORE2 > median(d$SCORE2), "high", "low"))
d$ZONE3 <- factor(ifelse(d$complexity > median(d$complexity), "hi", "lo"))
d2 <- FSSgam::case_study2
d2$Location <- factor(d2$Location); d2$Status <- factor(d2$Status)

scen <- list(
  cont_only        = d[, c("depth", "complexity", "SCORE2")],
  fact_only        = d[, c("ZONE", "ZONE2", "ZONE3")],
  mixed            = d[, c("depth", "complexity", "ZONE", "ZONE2")],
  one_factor       = d[, c("depth", "ZONE"), drop = FALSE],
  single_fact_col  = d[, "ZONE", drop = FALSE],
  single_cont_col  = d[, "depth", drop = FALSE],
  cs2_mixed        = na.omit(d2[, c("Distance", "Location", "Status")])
)

res <- list()
for (nm in names(scen)) {
  w <- NULL
  v <- withCallingHandlers(
    try(check_correlations(scen[[nm]]), silent = TRUE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
  )
  res[[nm]] <- list(value = v, warnings = w)
  nl <- withCallingHandlers(
    try(check_non_linear_correlations(scen[[nm]]), silent = TRUE),
    warning = function(x) invokeRestart("muffleWarning")
  )
  res[[paste0(nm, "__nonlinear")]] <- list(value = nl, warnings = NULL)
}
saveRDS(res, out.file)
cat("captured", length(res), "scenarios to", out.file, "\n")
