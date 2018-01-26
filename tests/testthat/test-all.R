rm(list=ls())
library(mediation)
library(testthat)
context("tests mediation")

require(parallel)
require(MASS)

accuracy <- 0.0001
accuracy2 <- 0.05

test_that("tests ivmediate on the jobs data", {
  # set random seed
  set.seed(12345)

  data(jobs)
 
  a <- lm(comply ~ treat + sex + age + marital + nonwhite + educ + income, data = jobs)
  b <- glm(job_dich ~ comply + treat + sex + age + marital + nonwhite + educ + income, data = jobs, family = binomial)
  c <- lm(depress2 ~ job_dich * (comply + treat) + sex + age + marital + nonwhite + educ + income, data = jobs)
  out <- ivmediate(a, b, c, sims = 50, boot = FALSE, enc = "treat", treat = "comply", mediator = "job_dich")
  
  expect_true("dc1.sims" %in% names(out))
  expect_equal(as.double(out$dc0.ci[1,1]), -0.07665396, tolerance = accuracy)
  expect_equal(as.double(out$dc1.ci["upper",1]), -0.00224326, tolerance = accuracy)
  expect_equal(as.double(out$dc1), -0.04293142, tolerance = accuracy)
  expect_equal(as.double(out$zc0), -0.04566892, tolerance = accuracy)
  expect_equal(out$dc1.sims[3], -0.06319376, tolerance = accuracy)
})  


test_that("tests mediator.design on the jobs data", {
  # set random seed
  set.seed(12345)

  data(jobs)
  foo.1 <- mediate.sed("depress2", "job_disc", "treat", jobs, SI=TRUE)
  x = summary(foo.1)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("d0.ci" %in% names(x))
  expect_equal(x$d0, -0.01557728, tolerance = accuracy)
  expect_equal(x$d1, -0.01706482, tolerance = accuracy)
  expect_equal(x$d0.ci[2], 0.01359242, tolerance = accuracy)
  expect_equal(x$d1.ci[1], -0.041206826, tolerance = accuracy)
  expect_true(x$design == "SED.NP.SI")

  foo.2 <- mediate.sed("depress2", "job_disc", "treat", jobs, SI=TRUE, boot=TRUE)
  x = summary(foo.2)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("d0.ci" %in% names(x))
  expect_equal(x$d0, -0.01557728, tolerance = accuracy)
  expect_equal(x$d1, -0.01706482, tolerance = accuracy)
  expect_equal(as.double(x$d0.ci[2]), 0.07442675, tolerance = accuracy)
  expect_equal(as.double(x$d1.ci[1]), -0.1759621, tolerance = accuracy)
  expect_true(x$design == "SED.NP.SI")
})  


test_that("tests mediate on the jobs data", {
  # set random seed
  set.seed(12345)
  
  data(jobs)
  ####################################################
  # Example 1: Linear Outcome and Mediator Models
  ####################################################
  b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
  c <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs)
  # Estimation via quasi-Bayesian approximation
  contcont <- mediate(b, c, sims=50, treat="treat", mediator="job_seek")
  out = summary(contcont)
  expect_true("tau.coef" %in% names(out))
  expect_equal(as.double(out$tau.coef), -0.05363842, tolerance = accuracy)
  expect_equal(out$d0.sims[3], -0.0323884, tolerance = accuracy)
  expect_equal(out$d.avg.sims[5], -0.03230468, tolerance = accuracy)
  expect_equal(out$n1.p, 0.32, tolerance = accuracy)
  expect_equal(as.double(out$model.y$coefficients['treat']), -0.0402647, tolerance = accuracy)
})  

