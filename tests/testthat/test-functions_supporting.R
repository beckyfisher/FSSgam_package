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

test_that("extract.mod.dat extracts AICc/BIC/r2 for a fitted gam model", {
  library(mgcv)
  data(case_study1)
  fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
    family = tw(), data = case_study1
  )

  out.r2 <- extract.mod.dat(fit, r2.type. = "r2")
  expect_equal(out.r2$AICc, MuMIn::AICc(fit), tolerance = 1e-6)
  expect_equal(out.r2$r2.vals, round(summary(fit)$r.sq, 5))

  out.dev <- extract.mod.dat(fit, r2.type. = "dev")
  expect_equal(out.dev$r2.vals, round(summary(fit)$dev.expl, 5))
})

test_that("extract.mod.dat returns an all-NA list for a try-error input", {
  out <- extract.mod.dat(structure("boom", class = "try-error"))
  expect_true(all(vapply(out, is.na, logical(1))))
  expect_named(out, c("AICc", "BIC", "r2.vals", "r2.vals.unique", "edf", "edf.less.1"))
})

test_that("build.inclusion.mat marks predictors present, regardless of interaction syntax", {
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
  mat <- build.inclusion.mat(included.vars, formula.list)

  # build.inclusion.mat does not set rownames -- it relies on row *position*
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

test_that("fit.mod.l updates a test.fit with a new formula on the supplied data", {
  library(mgcv)
  data(case_study1)
  base.fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
    family = tw(), data = case_study1
  )

  updated <- fit.mod.l(
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

test_that("fit.mod.l returns a try-error for a formula that cannot be fit", {
  library(mgcv)
  data(case_study1)
  base.fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
    family = tw(), data = case_study1
  )

  result <- fit.mod.l(
    formula.l = ~ s(not_a_column, k = 3, bs = "cr"),
    test.fit. = base.fit,
    use.dat = case_study1
  )
  expect_s3_class(result, "try-error")
})
