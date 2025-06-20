TODO list for the *beverstan* package
================
Yves Deville
2025-06-04

## BUG fixes

- [ ] Though not really bugs, fix the problems due to the re-definition
  of S3 generics such as `profLik`. This requires to change the
  dependences between packages.

## New features

- [ ] Implement the `show` method for the `"TVGEVBayes"` class.

- [ ] Implement `confInt` or `credInt` method for the `"TVGEVBayes"`
  class.

- [ ] Implement `residuals` (a.k.a; `resid`) method for the
  `"TVGEVBayes"` class to display the generalized residuals.

- [ ] What about some Information Criterion, e.g. WAIC (*Widely
  Applicable* or *Watababe Akaike*) Information Criterion?

- [ ] Implement the `simulate` method for the `"TVGEVBayes"` class. How
  to cope with the parameters? We could either plug the MAP parameter or
  sample the parameter from the its posterior, hence simulate from the
  predictive distribution.

- [ ] Allow the `predict` method for the `"TVGEVBayes"` class to compute
  the density or other “predict” things as done in the *revdbayes*
  package.

- [ ] Allow the `predict` method for the `"TVGEVBayes"` class to compute
  the distribution or density, … on *several* time ranges. The
  `autoplot` method could then plot e.g. the densities as is done in the
  *NSGEV* package.

- [ ] Implement the methods wanted in the **stan** packages such as
  `posterior_predict`.

- [ ] What about the implementation of random effects as in **INLA**:
  random walk, integregrated random walk, AR, …

## Known problems

- [ ] The initial values for the MCMC sampling may fail to match the
  data. These are obtained by applying a small perturbation on the
  parameters computed by ML and it may be the cas that the support of
  the GEV distributions does not cover all the observations, leading to
  an error. The solution would be to write a quite sophisticated
  function to get initial values.

- [ ] The method for the class may fail. The problem is that the
  quantiles of the predictive distribution are obtained by using on the
  which requires suitable lower and upper bounds for the root. These
  bounds are obtained by an heuristic which needs imporvements. The
  solution may be first to sample from the predictive distribution
  (which is easy) and simply take the maximum and the minimum of the
  simulated values as the bounds.

## Improvements

- [ ] The `autoplot` method for `TVGEVBayes` could have a `type`
  argument controlling the type of plot. The posssible values would
  be:  
  `fitted` (default) `RL`, `predict` and `residuals`. Depending of the
  value of `which` we would first call the corresponding method and the
  arguments passed through the dots `...` would be passed to the
  corresponding method. Each of these main methods would return an
  object with a corresponding class (e.g. `fitted.TVGVBayes`) inheriting
  from `"data.frame"` with some extra attributes to keep trace of the
  choices: `timeRange`, `level`, `what` (for what is predicted), … Note
  that this is some kind of refactoring.

- [ ] `autoplot` methods could provide a default title, such as
  `"posterior mean for the GEV expectation"`. Similarly, the predictive
  distribution could be described.

- [ ] Improve the management of the blocks. A duration of several years
  could be allowed. The beginning of the blocks could also be chosen
  differently if needed, for instance at the 1-st of August in order to
  allow a better management of winter extremes. This could be done by
  using a syntax like `blockBegin = "08-01"`.

- [ ] Define a class `"ibts"` for irregularly sampled block time series.
  Each observation would correspond to a time-block such as one year,
  and the blocks would be given in time order but with possible gaps. We
  could allow interval censored observations; the name of the variable
  would then no longer be a colname but be attached as an attribute to
  the object. So the columns could be named ̀“y”`,`“yL”`and`“yU”\`
  without having some complex management of suffixes and prefixes.

## Technical solutions

- [ ] Study the code of some packages like **brms** to see how to
  achieve gains in the use of **Stan**. Also instructions collected from
  the **rstantools** package documentation can help.
