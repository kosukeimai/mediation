context("ivmediate is reproducible under set.seed() with multicore")

data("jobs", package = "mediation")
a <- lm(comply ~ treat + sex + age + marital + nonwhite + educ + income,
        data = jobs)
b <- glm(job_dich ~ comply + treat + sex + age + marital + nonwhite + educ +
           income, data = jobs, family = binomial)
c <- lm(depress2 ~ job_dich * (comply + treat) + sex + age + marital +
          nonwhite + educ + income, data = jobs)

ivsims <- function(boot, multicore, mc.cores) {
  set.seed(2024)
  out <- suppressWarnings(
    ivmediate(a, b, c, sims = 10, boot = boot, enc = "treat",
              treat = "comply", mediator = "job_dich",
              multicore = multicore, mc.cores = mc.cores)
  )
  out[c("dc1.sims", "dc0.sims", "zc1.sims", "zc0.sims", "tauc.sims")]
}

for (boot in c(FALSE, TRUE)) {
  test_that(paste0("boot = ", boot, ": same draws serially and on 2 cores"), {
    # mclapply cannot fork on Windows; CRAN checks allow at most 2 cores.
    skip_on_os("windows")
    serial <- ivsims(boot, multicore = FALSE, mc.cores = 1L)
    expect_identical(ivsims(boot, multicore = TRUE, mc.cores = 2L), serial)
  })
}
