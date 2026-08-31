# Golden-master capture for the phase 3 restructure of the interaction-resolution
# code (specification section 5.3). Captures everything generate_model_set()
# returns that a restructure could plausibly disturb, for every argument
# combination, including ones the committed suite does not reach.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(mgcv)
out.file <- commandArgs(trailingOnly = TRUE)[1]

use.dat <- FSSgam::case_study1
use.dat$site <- as.factor(use.dat$site)
use.dat$ZONE2 <- factor(ifelse(use.dat$SCORE2 > median(use.dat$SCORE2), "high", "low"))
use.dat$ZONE3 <- factor(ifelse(use.dat$complexity > median(use.dat$complexity), "hi", "lo"))
# a deliberately near-duplicate factor, so the cov.cutoff screen actually bites
use.dat$ZONE.copy <- use.dat$ZONE
use.dat$ZONE4 <- factor(ifelse(use.dat$depth > median(use.dat$depth), "deep", "shallow"))
use.dat$ZONE5 <- factor(ifelse(use.dat$rugosity > median(use.dat$rugosity), "r1", "r2"))

cs3 <- FSSgam::case_study3
cs3$year <- factor(cs3$year); cs3$Sex <- factor(cs3$Sex); cs3$Species <- factor(cs3$Species)

test.fit <- gam(log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
                data = use.dat)
cs3.fit <- gam(GSI ~ s(month, k = 4, bs = "cc"), data = cs3)

base <- list(use.dat = use.dat, test.fit = test.fit,
             null.terms = "s(site,bs='re')", k = 3)

# a user-supplied correlation matrix over the three continuous predictors plus
# ZONE, saying everything is uncorrelated
cm.vars <- c("depth", "complexity", "SCORE2", "ZONE")
supplied.cm <- matrix(0, length(cm.vars), length(cm.vars),
                      dimnames = list(cm.vars, cm.vars))
diag(supplied.cm) <- 1

# and one over two factors, for the factor.factor.interactions screen
ff.vars <- c("ZONE", "ZONE.copy")
supplied.ff <- matrix(0, 2, 2, dimnames = list(ff.vars, ff.vars))
diag(supplied.ff) <- 1

scen <- list(
  cont_only = list(pred.vars.cont = c("complexity", "depth"), max.predictors = 2),
  fact_only = list(pred.vars.cont = NA, pred.vars.fact = c("ZONE", "ZONE2"),
                   max.predictors = 2),
  mixed = list(pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
               max.predictors = 2),
  ffi_true = list(pred.vars.cont = "depth", pred.vars.fact = c("ZONE", "ZONE2"),
                  factor.factor.interactions = TRUE, max.predictors = 2),
  ffi_true_3fact = list(pred.vars.cont = "depth",
                        pred.vars.fact = c("ZONE", "ZONE2", "ZONE3"),
                        factor.factor.interactions = TRUE, max.predictors = 3),
  ffi_true_screened = list(pred.vars.cont = "depth",
                           pred.vars.fact = c("ZONE", "ZONE.copy", "ZONE2"),
                           factor.factor.interactions = TRUE, max.predictors = 2),
  ffi_char = list(pred.vars.cont = "depth",
                  pred.vars.fact = c("ZONE", "ZONE2", "ZONE3"),
                  factor.factor.interactions = c("ZONE", "ZONE2"), max.predictors = 2),
  ffi_char_screened = list(pred.vars.cont = "depth",
                           pred.vars.fact = c("ZONE", "ZONE.copy", "ZONE2"),
                           factor.factor.interactions = c("ZONE", "ZONE.copy", "ZONE2"),
                           max.predictors = 2),
  ssi_true = list(pred.vars.cont = c("depth", "SCORE2", "complexity"),
                  smooth.smooth.interactions = TRUE, max.predictors = 2),
  ssi_true_max3 = list(pred.vars.cont = c("depth", "SCORE2", "complexity"),
                       smooth.smooth.interactions = TRUE, max.predictors = 3),
  ssi_true_screened = list(pred.vars.cont = c("depth", "SCORE2", "complexity"),
                           smooth.smooth.interactions = TRUE, cov.cutoff = 0.1,
                           max.predictors = 2),
  ssi_char = list(pred.vars.cont = c("depth", "SCORE2", "complexity"),
                  smooth.smooth.interactions = c("depth", "SCORE2"), max.predictors = 2),
  ssi_char_max3 = list(pred.vars.cont = c("depth", "SCORE2", "complexity"),
                       smooth.smooth.interactions = c("depth", "SCORE2", "complexity"),
                       max.predictors = 3),
  fsi_list = list(pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
                  factor.smooth.interactions = list(fact.vars = "ZONE",
                                                    cont.vars = "depth"),
                  max.predictors = 2),
  fsi_na = list(pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
                factor.smooth.interactions = NA, max.predictors = 2),
  fsi_char = list(pred.vars.cont = c("complexity", "depth"),
                  pred.vars.fact = c("ZONE", "ZONE2"),
                  factor.smooth.interactions = "ZONE", max.predictors = 2),
  cyclic = list(pred.vars.cont = c("month", "lunar.date"), cyclic.vars = "month",
                max.predictors = 2, k = 4,
                use.dat = cs3, test.fit = cs3.fit, null.terms = ""),
  linear_vars = list(pred.vars.cont = "depth", pred.vars.fact = "ZONE",
                     linear.vars = "complexity", max.predictors = 2),
  linear_vars_list_fsi = list(pred.vars.cont = "depth", pred.vars.fact = "ZONE",
                              linear.vars = "complexity",
                              factor.smooth.interactions = list(
                                fact.vars = "ZONE", cont.vars = "depth",
                                linear.vars = "complexity"),
                              max.predictors = 2),
  supplied_cor_matrix = list(pred.vars.cont = c("depth", "complexity", "SCORE2"),
                             pred.vars.fact = "ZONE", cor.matrix = supplied.cm,
                             smooth.smooth.interactions = TRUE, max.predictors = 2),
  supplied_cor_matrix_ffi = list(pred.vars.cont = "depth",
                                 pred.vars.fact = c("ZONE", "ZONE.copy"),
                                 factor.factor.interactions = TRUE,
                                 cor.matrix = supplied.ff, max.predictors = 2),
  non_linear_cor = list(pred.vars.cont = c("depth", "SCORE2", "complexity"),
                        smooth.smooth.interactions = TRUE,
                        non.linear.correlations = TRUE, max.predictors = 2),
  # Discriminating scenarios, added after measuring the baseline.
  #
  # The two smooth.smooth.interactions branches differ in arity: TRUE hard-codes
  # combinations of size 2, a character vector goes up to max.predictors. The
  # default cov.cutoff of 0.28 hides that here, because depth/complexity
  # correlate at 0.3336 and screening that pair out screens the three-way
  # combination containing it out as well. cov.cutoff = 0.4 admits all three
  # pairs, so a three-way te() becomes reachable through the character form.
  ssi_true_3way = list(pred.vars.cont = c("depth", "SCORE2", "complexity"),
                       smooth.smooth.interactions = TRUE, cov.cutoff = 0.4,
                       max.predictors = 3),
  ssi_char_3way = list(pred.vars.cont = c("depth", "SCORE2", "complexity"),
                       smooth.smooth.interactions = c("depth", "SCORE2", "complexity"),
                       cov.cutoff = 0.4, max.predictors = 3),

  # The two factor.factor.interactions branches differ in which half of the
  # correlation matrix they screen: the logical branch tests upper.tri and
  # lower.tri, the character branch only upper.tri. check_correlations() fits
  # multinom() separately in each direction, so the factor block is asymmetric
  # -- by at most 8.5e-4 over these columns. ZONE3/ZONE5 measure 0.5425 one way
  # and 0.5434 the other, so cov.cutoff = 0.543 falls between them and the two
  # branches must disagree about that pair.
  # A third factor is needed, uncorrelated with both. With only ZONE3 and ZONE5
  # the character branch skips the whole block through a separate guard,
  # length(which(factor.correlations < cov.cutoff)) > 1, which counts matrix
  # cells rather than pairs: one of the two off-diagonal cells is above the
  # cutoff, so the count is 1 and nothing is built, for a reason unrelated to
  # the screening difference under test.
  ffi_asym_screen = list(pred.vars.cont = "depth",
                         pred.vars.fact = c("ZONE", "ZONE3", "ZONE5"),
                         factor.factor.interactions = TRUE, cov.cutoff = 0.543,
                         max.predictors = 2),
  ffi_char_asym_screen = list(pred.vars.cont = "depth",
                              pred.vars.fact = c("ZONE", "ZONE3", "ZONE5"),
                              factor.factor.interactions = c("ZONE", "ZONE3", "ZONE5"),
                              cov.cutoff = 0.543, max.predictors = 2),
  # and the two-factor case kept, since it pins the cell-counting guard itself
  ffi_char_cellcount_guard = list(pred.vars.cont = "depth",
                                  pred.vars.fact = c("ZONE3", "ZONE5"),
                                  factor.factor.interactions = c("ZONE3", "ZONE5"),
                                  cov.cutoff = 0.543, max.predictors = 2),

  max_pred_1 = list(pred.vars.cont = "depth", pred.vars.fact = c("ZONE", "ZONE2"),
                    factor.factor.interactions = TRUE, max.predictors = 1),
  max_pred_1_ssi = list(pred.vars.cont = c("depth", "SCORE2"),
                        smooth.smooth.interactions = TRUE, max.predictors = 1),
  single_pred = list(pred.vars.cont = "depth", max.predictors = 1),
  single_fact = list(pred.vars.cont = NA, pred.vars.fact = "ZONE", max.predictors = 1)
)

deparse_one <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

res <- list()
for (nm in names(scen)) {
  # plain element replacement, not modifyList(): modifyList() recurses into
  # elements that are themselves lists, and a data.frame is a list, so it
  # merges use.dat column-by-column instead of replacing it
  args <- base
  args[names(scen[[nm]])] <- scen[[nm]]
  w <- NULL
  v <- withCallingHandlers(
    try(do.call(generate_model_set, args), silent = TRUE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
  )
  if (inherits(v, "try-error")) {
    res[[nm]] <- list(error = as.character(v), warnings = w)
  } else {
    res[[nm]] <- list(
      n.mods = v$n.mods,
      modnames = names(v$mod.formula),
      formulas = vapply(v$mod.formula, deparse_one, character(1)),
      predictor.correlations = v$predictor.correlations,
      used.data.names = names(v$used.data),
      used.data.dim = dim(v$used.data),
      included.vars = v$included.vars,
      warnings = w
    )
  }
}
saveRDS(res, out.file)
cat("captured", length(res), "scenarios\n")
for (nm in names(res)) {
  r <- res[[nm]]
  cat(sprintf("  %-24s %s  warnings=%d\n", nm,
      if (!is.null(r$error)) "ERROR" else paste0("n.mods=", r$n.mods),
      length(r$warnings)))
}
