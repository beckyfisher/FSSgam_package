# Tests for build_model_formulas(): how the resolved term names are turned into
# actual gam formulas. Everything here asserts on formula text rather than on
# fitted values.

# ---- cyclic.vars ------------------------------------------------------------

test_that("cyclic.vars substitutes bs='cc' into main-effect smooths", {
  fit <- fixture_cs3_cyclic()

  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = c("lunar.date", "month"),
    cyclic.vars = c("lunar.date", "month"),
    max.predictors = 1,
    k = 5
  )

  expect_equal(
    deparse_one(model.set$mod.formula[["lunar.date"]]),
    "~s(lunar.date, k = 5, bs = \"cc\")"
  )
  expect_equal(
    deparse_one(model.set$mod.formula[["month"]]),
    "~s(month, k = 5, bs = \"cc\")"
  )
})

test_that("cyclic.vars leaves non-cyclic continuous predictors on bs.arg", {
  fit <- fixture_cs3_cyclic()

  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = c("lunar.date", "gwt"),
    cyclic.vars = "lunar.date",
    max.predictors = 1,
    k = 5
  )

  expect_match(deparse_one(model.set$mod.formula[["lunar.date"]]), "bs = \"cc\"", fixed = TRUE)
  expect_match(deparse_one(model.set$mod.formula[["gwt"]]), "bs = \"cr\"", fixed = TRUE)
})

test_that("cyclic.vars substitutes bs='cc' into factor-smooth (by) terms", {
  fit <- fixture_cs3_cyclic()

  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = "lunar.date",
    pred.vars.fact = "Sex",
    cyclic.vars = "lunar.date",
    max.predictors = 2,
    k = 5
  )

  by.name <- grep(".by.", names(model.set$mod.formula), fixed = TRUE, value = TRUE)
  expect_length(by.name, 1)
  expect_match(
    deparse_one(model.set$mod.formula[[by.name]]),
    "s(lunar.date, by = Sex, k = 5, bs = \"cc\")",
    fixed = TRUE
  )
})

test_that("a te() of two cyclic predictors gets bs=c('cc','cc')", {
  fit <- fixture_cs3_cyclic()

  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = c("lunar.date", "month"),
    cyclic.vars = c("lunar.date", "month"),
    smooth.smooth.interactions = TRUE,
    max.predictors = 2,
    k = 5
  )

  expect_equal(
    deparse_one(model.set$mod.formula[["lunar.date.te.month"]]),
    "~te(lunar.date, month, k = 5, bs = c(\"cc\", \"cc\"))"
  )
})

test_that("a te() mixing a cyclic and a non-cyclic predictor gets bs=c('cc','cr')", {
  fit <- fixture_cs3_cyclic()

  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = c("lunar.date", "gwt"),
    cyclic.vars = "lunar.date",
    smooth.smooth.interactions = TRUE,
    max.predictors = 2,
    k = 5
  )

  # the te() marginal bases follow the order the predictors appear in the term
  expect_equal(
    deparse_one(model.set$mod.formula[["lunar.date.te.gwt"]]),
    "~te(lunar.date, gwt, k = 5, bs = c(\"cc\", \"cr\"))"
  )
})

test_that("cyclic model sets fit, and shrinkage is reflected in edf", {
  fit <- fixture_cs3_cyclic()

  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = c("lunar.date", "month"),
    pred.vars.fact = "Sex",
    cyclic.vars = c("lunar.date", "month"),
    max.predictors = 2,
    k = 5
  )
  out <- fit_quietly(model.set, parallel = FALSE)

  expect_length(out$failed.models, 0)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
  # bs='cc' uses shrinkage, so extract_mod_dat() resets any edf < 1 to 1 before
  # summing; every model must therefore report at least one parameter per term
  expect_true(all(out$mod.data.out$edf >= 1))
})

# ---- linear.vars ------------------------------------------------------------

test_that("linear.vars appear as bare linear terms, not smooths", {
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = "complexity",
    pred.vars.fact = NA,
    linear.vars = "SCORE2",
    max.predictors = 2
  )

  expect_equal(
    deparse_one(model.set$mod.formula[["SCORE2"]]),
    "~SCORE2 + s(site, bs = \"re\")"
  )
  expect_false(grepl("s(SCORE2", deparse_one(model.set$mod.formula[["SCORE2"]]), fixed = TRUE))
})

test_that("linear.vars interact with factors as a product term", {
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = "complexity",
    pred.vars.fact = "ZONE",
    linear.vars = "SCORE2",
    max.predictors = 2
  )

  expect_true("SCORE2.t.ZONE+ZONE" %in% names(model.set$mod.formula))
  expect_equal(
    deparse_one(model.set$mod.formula[["SCORE2.t.ZONE+ZONE"]]),
    "~SCORE2 * ZONE + ZONE + s(site, bs = \"re\")"
  )
})

test_that("a continuous predictor named in both pred.vars.cont and linear.vars is fitted linearly", {
  # setdiff(pred.vars.cont, linear.vars) decides which predictors get a smooth,
  # so naming one in both places makes it linear.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = NA,
    linear.vars = "depth",
    max.predictors = 1
  )

  expect_equal(
    deparse_one(model.set$mod.formula[["depth"]]),
    "~depth + s(site, bs = \"re\")"
  )
  expect_match(
    deparse_one(model.set$mod.formula[["complexity"]]), "s(complexity", fixed = TRUE
  )

  # This configuration is also the one that exercises enumerate_candidate_models()'s
  # deduplication: depth reaches use.mods twice, once from the pred.vars.cont
  # combn() and once from the linear.vars one, so without the unique() the set
  # carries two candidates both named "depth". Nothing downstream notices --
  # mod.formula[["depth"]] returns the first match either way -- but
  # fit_model_set() would fit and weight the same model twice.
  expect_equal(anyDuplicated(names(model.set$mod.formula)), 0L)
  expect_equal(
    anyDuplicated(vapply(model.set$mod.formula, deparse_one, character(1))), 0L
  )
})

test_that("a model set can be built from linear.vars alone", {
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = NA,
    pred.vars.fact = NA,
    linear.vars = c("complexity", "SCORE2"),
    max.predictors = 2
  )
  out <- fit_quietly(model.set, parallel = FALSE)

  expect_false(any(grepl("s(complexity", unlist(lapply(model.set$mod.formula, deparse_one)),
                          fixed = TRUE)))
  expect_length(out$failed.models, 0)
})

test_that("every linear.vars entry gets its own factor interaction term", {
  # Regression test: the pattern used to find .t. terms was built as
  # paste(linear.vars, ".t."), a vector whenever more than one linear predictor
  # was supplied. grep() then used only its first element, so the second and
  # later linear predictors' interaction terms were silently dropped from the
  # fitted formula -- "SCORE2.t.ZONE+ZONE" was named as an interaction but
  # fitted as the plain ZONE main-effect model.
  expect_no_warning(
    model.set <- fixture_cs1_model_set(
      pred.vars.cont = NA,
      pred.vars.fact = "ZONE",
      linear.vars = c("complexity", "SCORE2"),
      max.predictors = 2
    )
  )

  expect_equal(
    deparse_one(model.set$mod.formula[["ZONE+complexity.t.ZONE"]]),
    "~complexity * ZONE + ZONE + s(site, bs = \"re\")"
  )
  expect_equal(
    deparse_one(model.set$mod.formula[["SCORE2.t.ZONE+ZONE"]]),
    "~SCORE2 * ZONE + ZONE + s(site, bs = \"re\")"
  )
})

test_that("a list-form factor.smooth.interactions still fits its linear interaction", {
  # Regression test: the .t. terms are generated by resolve_factor_interactions()
  # from factor.smooth.interactions$linear.vars, which is validated against all
  # the predictors rather than against linear.vars. Deriving the match in
  # build_model_formulas() from linear.vars instead meant that where the two
  # disagreed the candidate was named as an interaction but fitted as the plain
  # factor main effect -- the same observable defect as the multiple-linear.vars
  # bug, reached by a different route.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = "ZONE",
    factor.smooth.interactions = list(
      fact.vars = "ZONE", cont.vars = "complexity", linear.vars = "depth"
    )
  )

  expect_true("ZONE+depth.t.ZONE" %in% names(model.set$mod.formula))
  expect_equal(
    deparse_one(model.set$mod.formula[["ZONE+depth.t.ZONE"]]),
    "~depth * ZONE + ZONE + s(site, bs = \"re\")"
  )
})

# ---- bs.arg -----------------------------------------------------------------

test_that("bs.arg controls the basis of every smooth term", {
  model.set <- fixture_cs1_model_set(bs.arg = "'ts'")
  formulas <- unlist(lapply(model.set$mod.formula, deparse_one))

  smooths <- formulas[grepl("s(complexity", formulas, fixed = TRUE) |
                             grepl("s(depth", formulas, fixed = TRUE)]
  expect_true(length(smooths) > 0)
  expect_true(all(grepl("bs = \"ts\"", smooths, fixed = TRUE)))
  expect_false(any(grepl("bs = \"cr\"", smooths, fixed = TRUE)))
})

test_that("a non-default bs.arg still produces a fittable model set", {
  model.set <- fixture_cs1_model_set(bs.arg = "'ts'")
  out <- fit_quietly(model.set, parallel = FALSE)

  expect_length(out$failed.models, 0)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
})

# ---- null.terms -------------------------------------------------------------

test_that("null.terms = '' gives an intercept-only null and no appended term", {
  model.set <- fixture_cs1_model_set(null.terms = "")

  expect_equal(deparse_one(model.set$mod.formula[["null"]]), "~1")
  expect_false(any(grepl("bs = \"re\"", unlist(lapply(model.set$mod.formula, deparse_one)),
                          fixed = TRUE)))
})

test_that("null.terms is appended to every candidate formula", {
  model.set <- fixture_cs1_model_set()
  formulas <- unlist(lapply(model.set$mod.formula, deparse_one))

  expect_true(all(grepl("s(site, bs = \"re\")", formulas, fixed = TRUE)))
})

# ---- gamm4 t2() -------------------------------------------------------------

test_that("a gamm4 test.fit produces t2() rather than te() smooth-smooth interactions", {
  # mgcv's te() is not available to gamm4, which needs a tensor-product basis
  # expressible as a sum of separable random effects; build_model_formulas()
  # therefore emits t2() when test.fit is of class gamm4.
  fit <- fixture_coral_ugamm()
  expect_s3_class(fit$test.fit, "gamm4")

  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    # av.wave and Depth correlate at 0.67 in coral_data, above the default
    # cov.cutoff, so their te() would be excluded before it could be built
    pred.vars.cont = c("Depth", "Survey"),
    smooth.smooth.interactions = TRUE,
    max.predictors = 2,
    k = 4
  )

  te.name <- grep(".te.", names(model.set$mod.formula), fixed = TRUE, value = TRUE)
  expect_length(te.name, 1)
  expect_match(deparse_one(model.set$mod.formula[[te.name]]), "^~t2\\(")
  expect_false(grepl("te(", deparse_one(model.set$mod.formula[[te.name]]), fixed = TRUE))
})

# ---- cyclic.vars are matched by name, not by substring ----------------------

test_that("a cyclic variable does not make a longer name cyclic too", {
  # A1. The basis used to be applied by rewriting assembled term strings, with
  # grep(cyclic.vars[r], term). That matches anywhere in the string, so
  # declaring "depth" cyclic also rewrote the term built for "depthx".
  use.dat <- fixture_cs1_data()
  use.dat$depthx <- use.dat$SCORE2

  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 4, bs = "cc") + s(site, bs = "re"),
    data = use.dat
  )
  model.set <- generate_model_set(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("depth", "depthx", "complexity"),
    cyclic.vars = "depth", null.terms = "s(site,bs='re')",
    max.predictors = 1, k = 4
  )

  f <- vapply(model.set$mod.formula, deparse_one, character(1))
  expect_match(f[["depth"]], "s(depth, k = 4, bs = \"cc\")", fixed = TRUE)
  expect_match(f[["depthx"]], "s(depthx, k = 4, bs = \"cr\")", fixed = TRUE)
  expect_match(f[["complexity"]], "s(complexity, k = 4, bs = \"cr\")", fixed = TRUE)
})

test_that("a full stop in a predictor name is not treated as a wildcard", {
  # The name was also used as a regular expression, so "." matched any
  # character and "a.b" made "axb" cyclic.
  use.dat <- fixture_cs1_data()
  use.dat$a.b <- use.dat$depth
  use.dat$axb <- use.dat$SCORE2

  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(a.b, k = 4, bs = "cc") + s(site, bs = "re"),
    data = use.dat
  )
  model.set <- generate_model_set(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("a.b", "axb"), cyclic.vars = "a.b",
    null.terms = "s(site,bs='re')", max.predictors = 1, k = 4
  )

  f <- vapply(model.set$mod.formula, deparse_one, character(1))
  expect_match(f[["a.b"]], "bs = \"cc\"", fixed = TRUE)
  expect_match(f[["axb"]], "bs = \"cr\"", fixed = TRUE)
})

test_that("a cyclic variable inside a te() pair gets its own basis", {
  # The old code needed a second block to repair the te() terms the first had
  # damaged. Both are gone; each marginal now carries its own basis.
  use.dat <- fixture_cs3_data()
  test.fit <- mgcv::gam(GSI ~ s(month, k = 4, bs = "cc"), data = use.dat)

  model.set <- generate_model_set(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("month", "lunar.date"), cyclic.vars = "month",
    smooth.smooth.interactions = TRUE, null.terms = "",
    max.predictors = 2, k = 4, cov.cutoff = 0.9
  )

  f <- vapply(model.set$mod.formula, deparse_one, character(1))
  te.name <- grep(".te.", names(f), fixed = TRUE, value = TRUE)
  expect_length(te.name, 1)
  # month is cyclic and lunar.date is not, in that order
  expect_match(f[[te.name]], "bs = c(\"cc\", \"cr\")", fixed = TRUE)
})

test_that("a three-way te() carries one basis per marginal", {
  # The repair block took only the first two variables of a te() term, so a
  # three-way te() was given two bs values. mgcv warns "bs wrong length and
  # ignored" and silently substitutes its own default, discarding both bs.arg
  # and any cyclic.vars for that term.
  fit <- fixture_cs1_gaussian()
  model.set <- generate_model_set(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("depth", "SCORE2", "complexity"),
    smooth.smooth.interactions = c("depth", "SCORE2", "complexity"),
    null.terms = "s(site,bs='re')", cov.cutoff = 0.4, max.predictors = 3, k = 3
  )

  f <- vapply(model.set$mod.formula, deparse_one, character(1))
  three.way <- grep("depth.te.SCORE2.te.complexity", names(f), fixed = TRUE, value = TRUE)
  expect_length(three.way, 1)
  expect_match(f[[three.way]], "bs = c(\"cr\", \"cr\", \"cr\")", fixed = TRUE)

  # and mgcv accepts it without the wrong-length warning
  expect_no_warning(
    mgcv::gam(stats::update(model.set$mod.formula[[three.way]],
                            log.Herbivore.biomass ~ .),
              data = model.set$used.data)
  )
})
