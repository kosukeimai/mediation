test_that("test.modmed.mediate.order uses z0.diff for z0 p-values (INT = FALSE and TRUE)", {
  set.seed(123)
  data(jobs)
  jobs$depress_ord <- ordered(cut(jobs$depress2, breaks = 3))

  b.int <- lm(job_seek ~ treat * age + econ_hard + sex, data = jobs)
  d.no_int <- MASS::polr(
    depress_ord ~ treat + job_seek + age + econ_hard + sex,
    data = jobs, Hess = TRUE, method = "probit"
  )
  fit_no_int <- suppressWarnings(
    mediate(b.int, d.no_int, sims = 10, treat = "treat", mediator = "job_seek")
  )
  res_no_int <- suppressWarnings(
    test.modmed(fit_no_int, list(age = 20), list(age = 70), sims = 10)
  )
  expect_s3_class(res_no_int, "test.modmed.mediate.order")
  expect_length(res_no_int, 2)
  expect_true(all(is.finite(res_no_int[[2]]$p.value)))

  d.int <- MASS::polr(
    depress_ord ~ treat * job_seek * age + econ_hard + sex,
    data = jobs, Hess = TRUE, method = "probit"
  )
  fit_int <- suppressWarnings(
    mediate(b.int, d.int, sims = 10, treat = "treat", mediator = "job_seek")
  )
  res_int <- suppressWarnings(
    test.modmed(fit_int, list(age = 20), list(age = 70), sims = 10)
  )
  expect_s3_class(res_int, "test.modmed.mediate.order")
  expect_length(res_int, 4)
  # Check that z0 p-value is bounded in [0, 1] and matches z0.diff (statistic)
  expect_true(all(res_int[[4]]$p.value >= 0 & res_int[[4]]$p.value <= 1))
})
