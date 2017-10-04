context("Compare new outputs to outputs from original `mediation` code")

# Only run tests on my local machine
if (isTRUE(unname(Sys.info()["user"])=="weihuang")) {

## ----------------------------------------------------------------------------
## Parameters
## ----------------------------------------------------------------------------

# Use the same `sims` as `gen_ref_objects`.
s <- 2000
tol <- 5/sqrt(s)

# Load reference objects
load("ref_objects.RData")

## ----------------------------------------------------------------------------
## Functions
## ----------------------------------------------------------------------------

reldiff <- function(x, y, tolerance = tol) {
  if (any(abs(x - y) < tol, !is.finite(x - y)))
    x - y
  else
    (x - y) / y
}

getelems <- function(x) {
  setdiff(
    names(x)[sapply(x, is.numeric)], 
    c(grep("^.*\\.sims|\\.ci|model*$", names(x), value = TRUE))
  )
}

## ----------------------------------------------------------------------------
## Framing example
## ----------------------------------------------------------------------------

set.seed(2014)
data("framing", package = "mediation")

## ----------------------------------------------------------------------------
## Standard case
## ----------------------------------------------------------------------------

cat("\nDoing Standard case...\n")
med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + treat + age + educ + gender + income,
               data = framing, family = binomial("probit"))
new.med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo",
                       boot = TRUE, sims = s, 
                       use_speed = TRUE)

elems <- getelems(med.out)
test_that("standard case returns equal output as reference", {
  expect_equal(
    lapply(new.med.out[elems], class), 
    lapply(med.out[elems], class)
  )
  expect_equal(
    unname(
      mapply(reldiff, unlist(new.med.out[elems]), unlist(med.out[elems]))), 
    rep(0, length(unlist(med.out[elems]))),
    tolerance = tol)
})

## ----------------------------------------------------------------------------
## Interaction between treatment and mediator
## ----------------------------------------------------------------------------

cat("\nDoing T-M interaction...\n")
med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo * treat + age + educ + gender + income,
               data = framing, family = binomial("probit"))
new.med.out.i <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo",
                         boot = TRUE, use_speed = TRUE,
                         sims = s)

elems <- getelems(med.out.i)
test_that("interaction between treatment and mediator returns 
  equal output as reference", {
  expect_equal(
    lapply(new.med.out.i[elems], class), 
    lapply(med.out.i[elems], class)
  )
  expect_equal(
    unname(
      mapply(reldiff, unlist(new.med.out.i[elems]), unlist(med.out.i[elems]))
    ), 
    rep(0, length(unlist(med.out.i[elems]))),
    tolerance = tol)
})


## ----------------------------------------------------------------------------
## Moderated mediators
## ----------------------------------------------------------------------------

cat("\nDoing Moderated mediators...\n")
med.fit <- lm(emo ~ treat * age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + treat * age + emo * age + educ + gender
               + income, data = framing, family = binomial("probit"))
new.med.age20 <- mediate(med.fit, out.fit, treat = "treat",
                         mediator = "emo", covariates = list(age = 20), 
                         boot = TRUE, use_speed = TRUE,
                         sims = s)
new.med.age60 <- mediate(med.fit, out.fit, treat = "treat",
                         mediator = "emo", covariates = list(age = 60), 
                         boot = TRUE, use_speed = TRUE,
                         sims = s)

elems <- getelems(med.age20)
test_that("moderated mediators returns equal output as reference", {
  expect_equal(
    lapply(new.med.age20[elems], class), 
    lapply(med.age20[elems], class)
  )
  expect_equal(
    unname(
      mapply(reldiff, unlist(new.med.age20[elems]), unlist(med.age20[elems]))
    ), 
    rep(0, length(unlist(med.age20[elems]))),
    tolerance = tol)
})
test_that("moderated mediators returns equal output as reference", {
  expect_equal(
    lapply(new.med.age60[elems], class), 
    lapply(med.age60[elems], class)
  )
  expect_equal(
    unname(
      mapply(reldiff, unlist(new.med.age60[elems]), unlist(med.age60[elems]))
    ), 
    rep(0, length(unlist(med.age60[elems]))),
    tolerance = tol)
})

## ----------------------------------------------------------------------------
## Non-binary treatments
## ----------------------------------------------------------------------------

cat("\nDoing Non-binary treatments...\n")
med.fit <- lm(emo ~ cond + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + cond + age + educ + gender + income,
               data = framing, family = binomial("probit"))
new.med23.out <- mediate(med.fit, out.fit, treat = "cond", mediator = "emo",
                        control.value = 2, treat.value = 3, 
                        boot = TRUE, use_speed = TRUE,
                        sims = s)

elems <- getelems(med23.out)
test_that("moderated mediators returns equal output as reference", {
  expect_equal(
    lapply(new.med23.out[elems], class), 
    lapply(med23.out[elems], class)
  )
  expect_equal(
    unname(
      mapply(reldiff, unlist(new.med23.out[elems]), unlist(med23.out[elems]))
    ), 
    rep(0, length(unlist(med23.out[elems]))),
    tolerance = tol)
})

## ----------------------------------------------------------------------------
## Binary outcome and ordered mediator (from example)
## ----------------------------------------------------------------------------

cat("\nDoing Binary outcome and ordered mediator...\n")

data("jobs", package = "mediation")
jobs$job_disc <- as.factor(jobs$job_disc)

med.fit <- polr(job_disc ~ treat + econ_hard + sex + age, data = jobs,
                method = "probit", Hess = TRUE)
out.fit <- glm(work1 ~ treat * job_disc + econ_hard + sex + age, data = jobs,
               family = binomial(link = "probit"))
new.med.polr <- mediate(med.fit, out.fit, 
                        treat = "treat", mediator = "job_disc",
                        boot = TRUE, sims = s, long = FALSE)

elems <- getelems(med.polr)
test_that("binary y/ordered mediator returns equal output as reference", {
  expect_equal(
    lapply(new.med.polr[elems], class), 
    lapply(med.polr[elems], class)
  )
  expect_equal(
    unname(
      mapply(reldiff, unlist(new.med.polr[elems]), unlist(med.polr[elems]))
    ), 
    rep(0, length(unlist(med.polr[elems]))),
    tolerance = tol)
})


## ----------------------------------------------------------------------------
## Multilevel example -- no bootstrap methods implemented
## ----------------------------------------------------------------------------

# data("student", package = "mediation")
# library(lme4)

# cat("\nDoing Multilevel case...\n")
# med.fit <- glmer(attachment ~ catholic + gender + income + pared + (1|SCH_ID),
#                  family = binomial(link = "logit"), data = student)
# out.fit <- glmer(fight ~ catholic * attachment +
#                          gender + income + pared + (1 + attachment|SCH_ID),
#                  family = binomial(link = "logit"), data = student)
# med.out.ml <- mediate(med.fit, out.fit, treat = "catholic", 
#                       mediator = "attachment", boot = FALSE, sims = s)

## ----------------------------------------------------------------------------
## Group-level treatment
## ----------------------------------------------------------------------------

# data("school", package = "mediation")

# cat("\nDoing Group-level treatment...\n")
# med.fit <- lm(smorale ~  free + catholic + coed, data = school)
# out.fit <- lmer(late ~ free + smorale + catholic + coed +
#                        gender + income + pared + (1|SCH_ID),
#                 data = student)
# med.out.ml2 <- mediate(med.fit, out.fit, treat = "free", mediator = "smorale",
#                        control.value = 3, treat.value = 4, boot = FALSE, 
#                        sims = s)

## ----------------------------------------------------------------------------
## IV mediate (TODO: implement boot)
## ----------------------------------------------------------------------------

# data("jobs", package = "mediation")

# cat("\nDoing IV mediate...\n")
# a <- lm(comply ~ treat + sex + age + marital + nonwhite + educ + income, 
#         data = jobs)
# b <- glm(job_dich ~ comply + treat + sex + age + marital + 
#           nonwhite + educ + income, data = jobs, family = binomial)
# c <- lm(depress2 ~ job_dich * (comply + treat) + sex + age + marital + 
#           nonwhite + educ + income, data = jobs)
# new.med.iv.out <- ivmediate(a, b, c, 
#                             enc = "treat", treat = "comply", 
#                             mediator = "job_dich",
#                             boot = TRUE, use_speed = TRUE,
#                             sims = s)

}