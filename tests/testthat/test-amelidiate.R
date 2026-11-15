test_that("amelidiate extracts n1, n0 from cis$n1.sims/n0.sims and passes z0 (not z1) to plot_process", {
  set.seed(123)
  data(jobs)
  d1 <- jobs[1:200, ]
  d2 <- jobs[201:400, ]
  datasets <- list(T1 = d1, T2 = d2)

  meds <- mediations(
    datasets,
    treatment = c("treat", "treat"),
    mediators = "job_seek",
    outcome = "depress2",
    covariates = "econ_hard + sex",
    families = c("gaussian", "gaussian"),
    interaction = TRUE,
    sims = 20
  )

  res <- amelidiate(meds)
  expect_false(is.null(res$n1))
  expect_false(is.null(res$n0))
  expect_equal(res$n1, mean(unlist(res$n1.sims)))
  expect_equal(res$n0, mean(unlist(res$n0.sims)))
  expect_equal(unname(res$coef.vec.0[2]), res$z0)
})
