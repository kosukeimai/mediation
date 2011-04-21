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
plot(out.3, treatment = "both")

model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital
                 + nonwhite + educ + income, data = jobs)
model.y <- gam(depress2 ~ treat + s(job_seek, bs = "cr")  + depress1 + econ_hard
                 + sex + age + occp + marital + nonwhite + educ + income, data = jobs)
out.4 <- mediate(model.m, model.y, sims = 1000, boot = TRUE, treat = "treat",
                  mediator = "job_seek")

model.y <- gam(depress2 ~ treat + s(job_seek, by = treat) + s(job_seek, by = control)
       + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income,
       data = jobs)
out.5 <- mediate(model.m, model.y, sims = 1000, boot = TRUE,
                    treat = "treat", mediator = "job_seek", control = "control")
summary(out.5)

model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital
              + nonwhite + educ + income, data = jobs)
model.y <- rq(depress2 ~ treat + job_seek + depress1 + econ_hard + sex
              + age + occp + marital + nonwhite + educ + income, 
              tau = 0.5, data = jobs)
out.6 <- mediate(model.m, model.y, sims = 1000, boot = TRUE, treat = "treat",
                  mediator = "job_seek")
summary(out.6)

model.y <- rq(depress2 ~ treat + job_seek + depress1 + econ_hard + sex
               + age + occp + marital + nonwhite + educ + income, 
               tau = 0.1, data = jobs)

model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age 
                + occp + marital + nonwhite + educ + income, data = jobs)
model.y <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age 
                + occp + marital + nonwhite + educ + income, 
                family = binomial(link = "probit"), data = jobs)
out.7 <- mediate(model.m, model.y, sims = 1000, treat = "treat",
                    mediator = "job_seek")
summary(out.7)

model.m <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age
                + occp + marital + nonwhite + educ + income, data = jobs, 
                family = binomial(link = "probit"))
model.y <- lm(depress2 ~ treat + job_dich + treat:job_dich + depress1
                + econ_hard + sex + age + occp + marital + nonwhite 
                + educ + income, data = jobs)
out.8 <- mediate(model.m, model.y, sims = 1000, treat = "treat",
                    mediator = "job_dich")
summary(out.8)

model.m <- polr(job_disc ~ treat + depress1 + econ_hard + sex + age 
                + occp + marital + nonwhite + educ + income, 
                data = jobs, method = "probit", Hess = TRUE)
model.y <- lm(depress2 ~ treat + job_disc + depress1 + econ_hard + sex + age 
                + occp + marital + nonwhite + educ + income, data = jobs)

out.9 <- mediate(model.m, model.y, sims = 1000, treat = "treat",
                     mediator = "job_disc")
summary(out.9)


model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp 
                + marital + nonwhite + educ + income, data = jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp
                + marital + nonwhite + educ + income, data = jobs)
med.cont <- mediate(model.m, model.y, sims=1000,  treat = "treat",
                       mediator = "job_seek")
                       
sens.cont <- medsens(med.cont, rho.by = 0.05)

summary(sens.cont)

plot(sens.cont, sens.par = "rho")

plot(sens.cont, sens.par = "R2", r.type = "total", sign.prod = "negative")

model.y <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age 
                + occp + marital + nonwhite + educ + income, 
                family = binomial(link = "probit"), data = jobs)
med.bout <- mediate(model.m, model.y, sims = 1000, treat = "treat",
                       mediator = "job_seek")

sens.bout <- medsens(med.bout, rho.by = 0.05, sims = 1000)

plot(sens.bout, sens.par = "rho")
plot(sens.bout, sens.par = "R2", r.type = "total", sign.prod = "positive")

plot(sens.bout, sens.par = "rho", pr.plot = TRUE)

model.m <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age 
                + occp + marital + nonwhite + educ + income, 
                data = jobs, family = binomial(link = "probit"))
model.y <- lm(depress2 ~ treat + job_dich+ depress1 + econ_hard + sex + age 
                + occp + marital + nonwhite + educ + income, data = jobs)
med.bmed <- mediate(model.m, model.y, sims = 1000, treat = "treat",
                       mediator = "job_dich")
sens.bmed <- medsens(med.bmed, rho.by = 0.05, sims = 1000)

plot(sens.bmed, sens.par = "rho")

