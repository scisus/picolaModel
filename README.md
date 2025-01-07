# picolaModel
<!-- badges: start -->
[![R-CMD-check](https://github.com/scisus/picolaModel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/scisus/picolaModel/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14597656.svg)](https://doi.org/10.5281/zenodo.14597656)
<!-- badges: end -->

This repository contains code to fit a multilevel Bayesian thermal time model
of pollen shed and cone receptivity for lodgepole pine (*Pinus contorta*
ssp. *latifolia*). The model uses [Stan](https://mc-stan.org) via the
['brms'](https://doi.org/10.32614/CRAN.package.brms) package.

For examples of making predictions using the model, see
[analysis/prediction_examples.md](analysis/prediction_examples.md).

For more insight into prior choices, see
[model/conceptualanalysisGDD.md](./model/conceptualanalysisGDD.md).

The model is fit to data from the
[`scisus/picolaDataFlowering`](https://github.com/scisus/picolaDataFlowering)
and [`scisus/picolaDataClimate`](https://github.com/scisus/picolaDataClimate)
packages.

This repository also includes diagnostic, output analysis and visualization
post-processing scripts that produce a variety of artifacts used by the
[`scisus/picolaDocThesis`](https://github.com/scisus/picolaDocThesis)
repository.

## Refitting and analysis

There are several possible ways to (re-)fit the model and/or run
analysis scripts. Scripts can be run interactively or from the command line.

```sh
# Run all preparation, fitting and analysis steps on fixed seed.
$ ./run_all.R --mode=fixed

# Same (--mode=fixed is the default)
$ ./run_all.R

# Run all preparation, fitting and analysis steps on a random seed.
$ ./run_all.R --mode=random

# Run all preparation and analysis steps but download prefit model.
# This might be preferable as the model is big and fitting is slow.
$ ./run_all.R --mode=download

# Just fit the model, don't run analysis scripts.
$ ./model/fit_model.R --mode=fixed

# Same (--mode=fixed is the default)
$ ./model/fit_model.R

# Just fit the model with a random seed.
$ ./model/fit_model.R --mode=random

# Just download the model object without refitting.
$ ./model/fit_model.R --mode=download

# Fit the model multiple times, randomly.
$ ./model/fit_model --mode=random --nruns=8

# Same, but perform diagnostics on each model fit.
$ ./model/fit_model --mode=random --nruns=8 --diagnostics=true
```
