## declare as Global Variables a number of variables that are used
## as intermediate variables in med.fun() and med.fun.ordered()

## These variables are created and assigned value in mediate()
## and are passed to med.fun() or med.fun.ordered() using
## environment(med.fun) <- environment() and
## environment(med.fun.ordered) <- environment()

if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("FamilyM", "cat.0", "cat.1", "control", "covariates", 
                           "isFactorM", "isFactorT", "isGlm.m", "isGlm.y", 
                           "isLm.m", "isLm.y", "isOrdered.m", "isRq.m", 
                           "isRq.y", "isSurvreg.m", "isSurvreg.y", "m",
                           "m.levels", "mediator", "model.m", "model.y", "n", 
                           "n.ycat", "outcome", "t.levels", "treat", 
                           "use_speed"))
}