test_that("wi computes Akaike weights that sum to 1", {
  w <- wi(c(100, 102, 105, 110))
  expect_equal(sum(w), 1, tolerance = 1e-8)
  # lowest AIC should get the highest weight
  expect_equal(which.max(w), 1L)
  expect_true(all(diff(w) <= 0))
})

test_that("wi gives equal weights when all AIC values are equal", {
  w <- wi(c(100, 100, 100))
  expect_equal(w, rep(1 / 3, 3))
})

test_that("wi propagates NA at its original position", {
  w <- wi(c(100, NA, 105))
  expect_true(is.na(w[2]))
  expect_false(anyNA(w[-2]))
  expect_equal(sum(w, na.rm = TRUE), 1, tolerance = 1e-8)
})

test_that("extract_mod_dat extracts AICc/BIC/r2 for a fitted gam model", {
  library(mgcv)
  data(case_study1)
  fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
    family = tw(), data = case_study1
  )

  out.r2 <- extract_mod_dat(fit, r2.type. = "r2")
  expect_equal(out.r2$AICc, MuMIn::AICc(fit), tolerance = 1e-6)
  expect_equal(out.r2$r2.vals, round(summary(fit)$r.sq, 5))

  out.dev <- extract_mod_dat(fit, r2.type. = "dev")
  expect_equal(out.dev$r2.vals, round(summary(fit)$dev.expl, 5))
})

test_that("extract_mod_dat returns an all-NA list for a try-error input", {
  out <- extract_mod_dat(structure("boom", class = "try-error"))
  expect_true(all(vapply(out, is.na, logical(1))))
  expect_named(out, c("AICc", "BIC", "r2.vals", "r2.vals.unique", "edf", "edf.less.1"))
})

test_that("build_inclusion_mat marks predictors present, regardless of interaction syntax", {
  included.vars <- c("a", "b", "c", "d")
  formula.list <- list(
    "null" = ~1,
    "a" = ~1,
    "a+b" = ~1,
    "a.by.b+b" = ~1,
    "a.I.b+a.I.b" = ~1,
    "a*b" = ~1,
    "a.t.b" = ~1,
    "a.te.b" = ~1
  )
  mat <- build_inclusion_mat(included.vars, formula.list)

  # build_inclusion_mat does not set rownames -- it relies on row *position*
  # matching the order of formula.list (fit_model_set() cbinds it onto
  # mod.data.out positionally), so index by position here too
  expect_equal(dim(mat), c(8, 4))
  expect_null(rownames(mat))
  expect_equal(colnames(mat), included.vars)
  expect_equal(unname(mat[1, ]), c(0, 0, 0, 0)) # "null"
  expect_equal(unname(mat[2, ]), c(1, 0, 0, 0)) # "a"
  for (i in 3:8) {
    expect_equal(unname(mat[i, ]), c(1, 1, 0, 0), info = names(formula.list)[i])
  }
})

test_that("fit_mod_l updates a test.fit with a new formula on the supplied data", {
  library(mgcv)
  data(case_study1)
  base.fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
    family = tw(), data = case_study1
  )

  updated <- fit_mod_l(
    formula.l = ~ s(complexity, k = 3, bs = "cr"),
    test.fit. = base.fit,
    use.dat = case_study1
  )

  expect_s3_class(updated, "gam")
  expect_match(as.character(formula(updated))[3], "complexity")
  # re-fitting the same spec independently should land at essentially the
  # same optimum (allow for Tweedie's iterative power-parameter estimation)
  manual <- update(base.fit, formula = ~ s(complexity, k = 3, bs = "cr"),
                    data = case_study1, family = family(base.fit))
  expect_equal(unname(coef(updated)), unname(coef(manual)), tolerance = 1e-4)
})

test_that("fit_mod_l keeps every refit's extended-family state independent, even when family is shared via a variable", {
  # Regression test for the conflict between GitHub issues #10 and #12.
  # family.opts[[2]] is the exact object gam() used to fit base.fit, so
  # fitting base.fit already mutates it in place (theta is stored in the
  # family object's own mutable environment). If refits of base.fit shared
  # that same object (or each other's), each one would retroactively change
  # an earlier refit's -- or base.fit's own -- already-estimated theta.
  #
  # family.opts is assigned to the global environment (via <<-) rather than
  # left local to this test_that() block -- see the matching comment in
  # test-fit_model_set.R for why that's required here regardless of this
  # fix (mgcv::gam() normalises away the environment a formula was
  # originally written in).
  library(mgcv)
  data(case_study1)
  use.dat <- case_study1
  use.dat$Herbivore.abundance <- round(use.dat$Herbivore.abundance)
  family.opts <<- list(gaussian(), nb())
  on.exit(rm("family.opts", envir = globalenv()), add = TRUE)
  base.fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
    family = family.opts[[2]], data = use.dat
  )
  theta.after.base.fit <- family.opts[[2]]$getTheta(TRUE)

  fit_a <- fit_mod_l(formula.l = ~ s(complexity, k = 3, bs = "cr"), test.fit. = base.fit, use.dat = use.dat)
  theta.a <- fit_a$family$getTheta(TRUE)

  fit_b <- fit_mod_l(formula.l = ~ s(depth, k = 3, bs = "cr") + s(complexity, k = 3, bs = "cr"),
                      test.fit. = base.fit, use.dat = use.dat)

  expect_equal(fit_a$family$getTheta(TRUE), theta.a, tolerance = 1e-10)
  expect_equal(family.opts[[2]]$getTheta(TRUE), theta.after.base.fit, tolerance = 1e-10)
  expect_false(identical(environment(fit_a$family$getTheta), environment(fit_b$family$getTheta)))
})

test_that("fit_mod_l returns a try-error for a formula that cannot be fit", {
  library(mgcv)
  data(case_study1)
  base.fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
    family = tw(), data = case_study1
  )

  result <- fit_mod_l(
    formula.l = ~ s(not_a_column, k = 3, bs = "cr"),
    test.fit. = base.fit,
    use.dat = case_study1
  )
  expect_s3_class(result, "try-error")
})

# ---- extract_mod_dat across model classes -----------------------------------

test_that("extract_mod_dat computes r2.lm.est for a gam from observed against fitted", {
  fit <- mgcv::gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
    family = tw(), data = FSSgam::case_study1
  )

  out <- extract_mod_dat(fit, r2.type. = "r2.lm.est")
  expected <- round(summary(stats::lm(fit$y ~ stats::predict(fit)))$r.sq, 5)

  expect_equal(out$r2.vals, expected)
  # r2.lm.est is computed on the link scale, so it is not the same quantity as
  # either summary()$r.sq or the deviance explained
  expect_false(isTRUE(all.equal(out$r2.vals, round(summary(fit)$r.sq, 5))))
})

test_that("extract_mod_dat sums edf plus parametric coefficients", {
  fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + ZONE,
    data = FSSgam::case_study1
  )

  out <- extract_mod_dat(fit, r2.type. = "r2")
  edf.m <- summary(fit)$edf
  edf.m[edf.m < 1] <- 1
  expect_equal(out$edf, round(sum(c(edf.m, length(summary(fit)$p.coeff))), 2))
  expect_equal(out$edf.less.1, length(which(summary(fit)$edf < 0.25)))
})

test_that("extract_mod_dat resets shrunk edf to 1 before summing", {
  # bs='cc' shrinks unsupported smooth terms towards zero edf; leaving those
  # below 1 would understate the model's parameter count in the output table
  use.dat <- FSSgam::case_study3
  fit <- mgcv::gam(
    GSI ~ s(lunar.date, k = 5, bs = "cc"), family = stats::Gamma(), data = use.dat
  )

  out <- extract_mod_dat(fit, r2.type. = "r2")
  expect_true(out$edf >= 1 + length(summary(fit)$p.coeff))
})

test_that("extract_mod_dat handles a gamm4 fit", {
  fit <- fixture_coral_ugamm()$test.fit
  expect_s3_class(fit, "gamm4")

  out.r2 <- extract_mod_dat(fit, r2.type. = "r2")
  expect_equal(out.r2$r2.vals, round(summary(fit$gam)$r.sq, 5))
  expect_equal(out.r2$AICc, MuMIn::AICc(fit), tolerance = 1e-6)
  expect_equal(out.r2$edf, round(sum(c(
    pmax(fit$gam$edf, 1), length(fit$gam$p.coeff)
  )), 2))

  # gamm4 does not report a deviance, so 'dev' is documented to be unavailable
  expect_true(is.na(extract_mod_dat(fit, r2.type. = "dev")$r2.vals))
})

test_that("extract_mod_dat computes r2.lm.est for a binomial gamm4 from the proportion", {
  # the binomial branch reconstructs the observed proportion from the
  # cbind(successes, failures) response before regressing it on the fitted values
  fit <- fixture_coral_ugamm()$test.fit

  out <- extract_mod_dat(fit, r2.type. = "r2.lm.est")
  y.dat <- attributes(fit$mer)$frame$y
  y <- y.dat[, 1] / (y.dat[, 1] + y.dat[, 2])
  x <- stats::predict(fit[[1]], re.form = NA, type = "response")
  expect_equal(out$r2.vals, round(summary(stats::lm(y ~ x))$r.sq, 5))
})

test_that("extract_mod_dat handles a gamm (nlme) fit", {
  fit <- fixture_coral_gamm()$test.fit
  expect_s3_class(fit, "gamm")

  out <- extract_mod_dat(fit, r2.type. = "r2")
  expect_equal(out$r2.vals, round(summary(fit$gam)$r.sq, 5))
  expect_equal(out$AICc, MuMIn::AICc(fit), tolerance = 1e-6)
  expect_true(is.finite(out$BIC))

  # Documented current behaviour, not an endorsement: extract_mod_dat() has no
  # gamm branch for the package default r2.type = "r2.lm.est", nor for "dev",
  # so both silently report NA for a uGamm(lme4 = FALSE) fit. Reported upstream
  # as its own issue; this test exists so that fixing it is a visible change.
  expect_true(is.na(extract_mod_dat(fit, r2.type. = "r2.lm.est")$r2.vals))
  expect_true(is.na(extract_mod_dat(fit, r2.type. = "dev")$r2.vals))
})

# ---- fit_mod_l --------------------------------------------------------------

test_that("fit_mod_l refits a gamm4 test.fit", {
  fit <- fixture_coral_ugamm()

  updated <- fit_mod_l(
    formula.l = cbind(successes, failures) ~ s(av.wave, k = 4, bs = "cr"),
    test.fit. = fit$test.fit,
    use.dat = fit$use.dat
  )

  expect_s3_class(updated, "gamm4")
  expect_match(deparse1(stats::formula(updated$gam)), "av.wave", fixed = TRUE)
})

test_that("fit_mod_l accepts an explicitly supplied family", {
  # fit_model_set() resolves the family once per candidate on the calling
  # process and passes it in, rather than relying on fit_mod_l()'s own default
  fit <- fixture_cs1_nb()
  fam <- resolve_candidate_family(fit$test.fit)

  updated <- fit_mod_l(
    formula.l = ~ s(complexity, k = 3, bs = "cr"),
    test.fit. = fit$test.fit, use.dat = fit$use.dat, family. = fam
  )

  expect_s3_class(updated, "gam")
  expect_identical(environment(updated$family$getTheta), environment(fam$getTheta))
})

test_that("fit_mod_l with family. = NULL leaves update() to its own behaviour", {
  # NULL is what resolve_candidate_family() returns for a test.fit whose call
  # carries no family=, i.e. the default gaussian()
  fit <- fixture_cs1_gaussian()

  updated <- fit_mod_l(
    formula.l = ~ s(complexity, k = 3, bs = "cr") + s(site, bs = "re"),
    test.fit. = fit$test.fit, use.dat = fit$use.dat, family. = NULL
  )

  expect_s3_class(updated, "gam")
  expect_equal(stats::family(updated)$family, "gaussian")
})

# fit_mod_l()'s dsm branch (length(grep("dsm", class(test.fit.))) > 0) is left
# uncovered deliberately. Exercising it needs a fitted object of class "dsm"
# from package dsm, which is neither an Import nor a Suggest here; adding a
# dependency on it -- and on the Distance/mrds stack it pulls in -- to cover two
# lines is not justified. The branch is a plain update() call with no data=
# argument, structurally identical to the covered branch.

test_that("extract_mod_dat computes r2.lm.est for a non-binomial gamm4 fit", {
  # the binomial branch reconstructs a proportion from cbind(successes,
  # failures); every other family takes the model frame's response directly
  use.dat <- fixture_coral_data()
  use.dat$prop <- use.dat$allcoral / use.dat$totalpoints
  fit <- MuMIn::uGamm(
    prop ~ s(Depth, k = 4, bs = "cr"), random = ~ (1 | Site),
    data = use.dat, lme4 = TRUE
  )
  expect_equal(stats::family(fit$mer)[[1]], "gaussian")

  out <- extract_mod_dat(fit, r2.type. = "r2.lm.est")
  expected <- round(summary(stats::lm(
    attributes(fit$mer)$frame$y ~
      stats::predict(fit[[1]], re.form = NA, type = "response")
  ))$r.sq, 5)

  expect_equal(out$r2.vals, expected)
  expect_false(is.na(out$r2.vals))
})
