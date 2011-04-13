################################################################################
# This file contains the R code used in the mediation vignette.
################################################################################
library("mediation")
data("jobs")
set.seed(4646)

model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age
                  + occp + marital + nonwhite + educ + income, data = jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age
                  + occp + marital + nonwhite + educ + income, data = jobs)

out.1 <- mediate(model.m, model.y, sims = 1000, boot = TRUE, treat = "treat",
                    mediator = "job_seek")
out.2 <- mediate(model.m, model.y, sims = 1000, treat = "treat",
                    mediator = "job_seek")

summary(out.2)

plot(out.2)

model.y <- lm(depress2 ~ treat + job_seek + treat:job_seek + depress1 + econ_hard
               + sex + age + occp + marital + nonwhite + educ + income, data = jobs)
out.3 <- mediate(model.m, model.y, sims = 1000, treat = "treat", 
                    mediator = "job_seek")
summary(out.3)

out.4 <- mediate(model.m, model.y, sims = 1000, boot = TRUE, treat = "treat",
                  mediator = "job_seek")

model.y <- gam(depress2 ~ treat + s(job_seek, by = treat) + s(job_seek, by = control)
       + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income,
       data = jobs)
       
out.5 <- mediate(model.m, model.y, sims = 1000, boot = TRUE,
                    treat = "treat", mediator = "job_seek", control = "control")
summary(out.5)

