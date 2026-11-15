test_that("plot.mediations plots the selected element in `which`, not x[[1]]", {
  set.seed(123)
  data(jobs)
  d1 <- jobs[1:200, ]
  datasets <- list(treat = d1)

  meds <- mediations(
    datasets,
    treatment = "treat",
    mediators = c("job_seek", "depress1"),
    outcome = "depress2",
    covariates = "econ_hard + sex",
    families = c("gaussian", "gaussian"),
    sims = 10
  )

  # Spy on plot.mediate calls or check plotted data via recordPlot / mock
  plotted <- character()
  with_mocked_bindings(
    plot.mediate = function(x, xlab = NULL, ...) {
      plotted <<- c(plotted, x$mediator)
    },
    {
      plot.mediations(meds, which = "depress2.treat.depress1", ask = FALSE)
    },
    .package = "mediation"
  )

  expect_equal(plotted, "depress1")
})
