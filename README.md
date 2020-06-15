# Leave-One-Covariate-Out Conditional Randomization Test (LOCO CRT)

The LOCO CRT is a computationally efficient version of the conditional randomization test that sacrifices little or no power compared to the CRT.

This repository contains code to reproduce the numerical simulations from the paper:

E. Katsevich and A. Ramdas, The leave-one-covariate-out conditional randomization test, 2020 (arXiv).

The script reproduce_all.R carries out our numerical experiments sequentially and plots the results. We do not recommend executing all of it directly, since our simulations are computationally intensive. We carried out these computations on the [Bridges supercomputer](https://www.psc.edu/bridges) through the NSF [XSEDE](https://www.xsede.org/) project. 

The LOCO CRT methodology is implemented in CRT.R, along with other CRT variants tested. The script run_one_experiment.R contains the main simulation code, with parameters specified by the two files starting with "input_file". 

## Dependencies

The following R (version 3.6.3) packages are used:

* glmnet 3.0
* reshape2 1.4.3
* MASS 7.3-51.4
* knockoff 0.3.2
* R.utils 2.9.0
* tidyverse 1.2.1