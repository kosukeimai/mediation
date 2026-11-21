test_that("additive lm effects agree draw by draw with the coefficient contrasts", {
  set.seed(81)
  d <- data.frame(t = rnorm(90), x = rnorm(90), g = factor(rep(1:3, 30)))
  d$w <- rep(c(0, 1, 3), length.out = nrow(d))
  d$m <- 3 + 0.8 * d$t - d$x + rnorm(nrow(d))
  d$y <- 1e9 + 2 * d$m - 0.3 * d$t + d$x + rnorm(nrow(d))
  d$w[1] <- 1
  m <- lm(m ~ x + t + g, data = d, weights = w)
  y <- lm(y ~ g + m + t + x, data = d, weights = w)
  B <- 60L
  set.seed(43)
  fit <- mediate(m, y, treat = "t", mediator = "m", control.value = -0.5,
                 treat.value = 0.25, sims = B)
  next_rng <- .Random.seed

  set.seed(43)
  alpha <- mvtnorm::rmvnorm(B, coef(m), vcov(m))
  beta <- mvtnorm::rmvnorm(B, coef(y), vcov(y))
  invisible(rnorm(B * nrow(d), sd = summary(m)$sigma))
  expect_identical(.Random.seed, next_rng)
  indirect <- 0.75 * alpha[, "t"] * beta[, "m"]
  direct <- 0.75 * beta[, "t"]
  total <- indirect + direct
  for (field in c("d0.sims", "d1.sims"))
    expect_equal(as.numeric(fit[[field]]), indirect, tolerance = 1e-12)
  for (field in c("z0.sims", "z1.sims"))
    expect_equal(as.numeric(fit[[field]]), direct, tolerance = 1e-12)
  expect_equal(as.numeric(fit$tau.sims), total, tolerance = 1e-12)
  expect_equal(as.numeric(fit$n0.sims), indirect / total, tolerance = 1e-12)
  expect_equal(unname(fit$d0.ci), unname(quantile(indirect, c(.025, .975))), tolerance = 1e-12)
  expect_equal(fit$d0.p, 2 * min(mean(indirect < 0), mean(indirect > 0)))
})

test_that("interactions, transformed variables, offsets and changed contrasts use predictions", {
  set.seed(20)
  d <- data.frame(t = rnorm(40), m = rnorm(40), y = rnorm(40),
                  g = factor(rep(letters[1:4], each = 10)))
  options_before <- options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(options_before))
  eligible <- function(fit) additive_lm_contrast(fit, model.frame(fit), c("t", "m"))
  expect_true(eligible(lm(y ~ g + t + m, d)))
  expect_false(eligible(lm(y ~ t * m, d)))
  expect_false(eligible(lm(y ~ t + m + I(m^2), d)))
  expect_false(eligible(lm(y ~ t + m + I(t^2), d)))
  expect_false(eligible(lm(y ~ t + m + offset(t), d)))
  expect_false(eligible(lm(y ~ t + m, d, offset = t)))
  expect_false(eligible(lm(y ~ g + t + m, d, contrasts = list(g = "contr.sum"))))
  expect_false(eligible(lm(y ~ t + m, transform(d, t = factor(t > 0)))))
  expect_false(eligible(lm(y ~ t + m, transform(d, m = factor(m > 0)))))

  # Check the effect of t:m against the fitted interaction coefficient.
  m <- lm(m ~ t, d)
  y <- lm(y ~ t * m, d)
  B <- 40L
  set.seed(8)
  fit <- mediate(m, y, treat = "t", mediator = "m", sims = B)
  set.seed(8)
  a <- mvtnorm::rmvnorm(B, coef(m), vcov(m))[, "t"]
  b <- mvtnorm::rmvnorm(B, coef(y), vcov(y))
  expect_equal(as.numeric(fit$d0.sims), a * b[, "m"], tolerance = 1e-12)
  expect_equal(as.numeric(fit$d1.sims), a * (b[, "m"] + b[, "t:m"]), tolerance = 1e-12)
})

test_that("factor indicators cannot shadow treatment or mediator coefficients", {
  set.seed(8)
  d <- data.frame(a = factor(rep(c("a", "b"), 45)), ab = rnorm(90), m = rnorm(90))
  d$m <- d$m + .7 * d$ab
  d$y <- 3 * d$ab + 8 * (d$a == "b") + 2 * d$m + rnorm(90)
  m <- lm(m ~ a + ab, d)
  y <- lm(y ~ a + ab + m, d)
  expect_false(additive_lm_contrast(m, model.frame(m), "ab"))
  expect_false(additive_lm_contrast(y, model.frame(y), c("ab", "m")))
  set.seed(5)
  fit <- mediate(m, y, treat = "ab", mediator = "m", sims = 30)
  set.seed(5)
  alpha <- mvtnorm::rmvnorm(30, coef(m), vcov(m))
  beta <- mvtnorm::rmvnorm(30, coef(y), vcov(y))
  # Each third column is the numeric treatment; the second column is the factor indicator.
  expect_equal(as.numeric(fit$z0.sims), beta[, 3], tolerance = 1e-12)
  expect_equal(as.numeric(fit$d0.sims), alpha[, 3] * beta[, "m"], tolerance = 1e-12)
})
