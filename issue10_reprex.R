# Reprex for https://github.com/beckyfisher/FSSgam/issues/10
# "Models fail to fit when parallel = TRUE"
#
# Reported behaviour: when `family` is passed to test.fit as an *element of
# a vector/list* (rather than written out literally, e.g. family = tw()),
# every candidate model fails to fit once parallel = TRUE -- but the same
# call set succeeds with parallel = FALSE.
#
# Root cause: fit_mod_l() (R/functions_supporting.R) refits each candidate
# via update(test.fit., formula = formula.l, data = use.dat), with no
# explicit family= override. update() therefore re-evaluates the *original*
# unevaluated family expression captured in test.fit$call (here,
# `family.vec[[2]]`). On the master process that expression resolves fine
# (family.vec is in scope), but under parallel = TRUE each candidate is
# fit on a fresh doSNOW worker process that never had family.vec exported
# to it, so the re-evaluation fails with "object 'family.vec' not found".
#
# NOTE: this exact line was *intentionally* changed on 2026-06-22 (commit
# 29daf70, fixing GitHub issue #12) to stop passing
# family = stats::family(test.fit.) into update(). That explicit-family
# pattern was itself added on 2022-06-29 (commit a0a901d, "fix for #10")
# specifically to solve *this* issue -- substituting an already-evaluated
# family object into the call sidesteps any re-evaluation of family.vec on
# the worker. Removing it (to stop mgcv's extended families warm-starting
# theta across refits, see issue #12) reopened this issue. The two bugs
# cannot both be fixed by simply toggling that one override on or off --
# see the discussion on #10/#12 for why.
#
# Run with `devtools::load_all()` rather than `library(FSSgam)`: an
# installed copy may be stale relative to the current source.

devtools::load_all(".", quiet = TRUE)
library(mgcv)

set.seed(1)

data(case_study1)
use.dat <- case_study1
use.dat$site <- as.factor(use.dat$site)

# The key ingredient: family stored in a vector/list and indexed into,
# instead of being written directly in the gam() call.
family.vec <- list(gaussian(), tw())
test.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
                 family = family.vec[[2]], data = use.dat)

model.set <- generate_model_set(
  use.dat = use.dat,
  test.fit = test.fit,
  pred.vars.cont = c("complexity", "depth"),
  pred.vars.fact = "ZONE",
  null.terms = "s(site,bs='re')",
  max.predictors = 2,
  k = 3
)

cat("\n=== parallel = FALSE ===\n")
fit.seq <- fit_model_set(model.set, parallel = FALSE)
cat("Models that failed to fit (sequential):", length(fit.seq$failed.models), "\n")
cat("Models that succeeded (sequential):", length(fit.seq$success.models), "\n")

cat("\n=== parallel = TRUE ===\n")
fit.par <- fit_model_set(model.set, parallel = TRUE, n.cores = 2)
cat("Models that failed to fit (parallel):", length(fit.par$failed.models), "\n")
cat("Models that succeeded (parallel):", length(fit.par$success.models), "\n")
