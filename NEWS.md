# mediation 4.5.0

A minor update that introduces a number of additional features.

## New features
- A new function `mediate_tsls()` that conducts mediation analysis using an instrumental variable estimator
- In `mediate()`, the `boot = TRUE` option now makes use of the `boot()` function in the `boot` package to implement nonparametric Bayesian bootstrap. Users may also pass arguments to the `boot()` function in the `boot` package e.g. `parallel` or `ncpus`.
- For `mediate()`, add a `use_speed` option to re-fit `lm` and `glm` models using functions from the `speedglm` package when nonparametric bootstrap is used.

## Bug fixes
- Fixed a bug that may cause the package to fail to load when the `system.file()` function from `base` is written over by other packages e.g. `pkgload`.