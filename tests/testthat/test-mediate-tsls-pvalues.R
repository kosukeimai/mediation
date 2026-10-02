context("Two-sided analytic p-values for mediate_tsls")

tsls_fixture <- function(direct = 0.2, mediated = 0.3, treatment = 1,
                         outcome_sign = 1) {
  n <- 64
  x <- seq_len(n)
  basis <- qr.Q(qr(cbind(1, x, x^2, x^3, x^4)))[, 2:5] * sqrt(n)
  dat <- data.frame(t = basis[, 1], z = basis[, 2])
  dat$m <- treatment * dat$t + 2 * dat$z + 0.5 * basis[, 3]
  dat$y <- outcome_sign * (direct * dat$t + mediated * dat$m + basis[, 4])
  model.m <- do.call(lm, list(formula = m ~ t + z, data = dat))
  model.y <- do.call(lm, list(formula = y ~ t + m, data = dat))
  mediate_tsls(model.m, model.y, treat = "t", boot = FALSE)
}

expect_two_sided_tsls <- function(out) {
  for (effect in c("d1", "z0", "tau", "n0")) {
    estimate <- if (effect == "tau") out$tau.coef else out[[effect]]
    standard_error <- out[[paste0(effect, ".se")]]
    p_value <- out[[paste0(effect, ".p")]]
    interval <- out[[paste0(effect, ".ci")]]
    expected <- 2 * pnorm(abs(estimate / standard_error), lower.tail = FALSE)
    expect_equal(unname(p_value), unname(expected), tolerance = 1e-12)
    expect_true(is.finite(p_value) && p_value >= 0 && p_value <= 1)
    expect_identical(unname(p_value < 1 - out$conf.level),
                     unname(interval[1] > 0 || interval[2] < 0))
  }
  expect_equal(out$d0.p, out$d1.p)
  expect_equal(out$z1.p, out$z0.p)
}

test_that("all four analytic p-values use both normal tails", {
  expect_two_sided_tsls(tsls_fixture())
})

test_that("analytic p-values agree with confidence intervals near 5 percent", {
  preliminary <- tsls_fixture()
  out <- tsls_fixture(direct = 1.8 * preliminary$z0.se)
  expect_equal(unname(out$z0 / out$z0.se), 1.8, tolerance = 1e-12)
  expect_equal(unname(out$z0.p), 2 * pnorm(-1.8), tolerance = 1e-12)
  expect_gt(out$z0.p, 0.05)
  expect_true(out$z0.ci[1] < 0 && out$z0.ci[2] > 0)
  expect_two_sided_tsls(out)
})

test_that("reversing outcome signs preserves two-sided p-values", {
  positive <- tsls_fixture()
  negative <- tsls_fixture(outcome_sign = -1)
  expect_two_sided_tsls(negative)
  for (effect in c("d1", "z0", "tau", "n0")) {
    expect_equal(positive[[paste0(effect, ".p")]],
                 negative[[paste0(effect, ".p")]], tolerance = 1e-12)
  }
})

test_that("zero effects have analytic p-values of one", {
  no_direct <- tsls_fixture(direct = 0)
  expect_equal(unname(no_direct$z0), 0, tolerance = 1e-12)
  expect_equal(unname(no_direct$z0.p), 1, tolerance = 1e-12)
  no_mediation <- tsls_fixture(treatment = 0)
  expect_equal(unname(no_mediation$d1), 0, tolerance = 1e-12)
  expect_equal(unname(no_mediation$d1.p), 1, tolerance = 1e-12)
  expect_equal(unname(no_mediation$n0.p), 1, tolerance = 1e-12)
  expect_two_sided_tsls(no_direct)
  expect_two_sided_tsls(no_mediation)
})
