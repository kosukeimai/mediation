fit_speedglm <- function(x) {
  speedglm::speedglm(formula = x$formula,
                     data = x$data,
                     family = eval(x$family),
                     weights = x$weights)
}
