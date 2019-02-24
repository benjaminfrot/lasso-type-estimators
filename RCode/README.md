# *Note* This code is deprecated. See: https://github.com/benjaminfrot/lrpsadmm .


This directory contains an R implementation of "Efficient Latent Variable Graphical Model Selection
via Split Bregman Method" (http://arxiv.org/abs/1110.3076), where the Alternative Direction Method of Multipliers (ADMM) is used to compute the estimator suggested in "Latent variable graphical model selection via convex optimization" (https://projecteuclid.org/euclid.aos/1351602527).

See test_split_bregman_lrps.R for how to use the code. Note that you might need to change the "source" statements at the top depending on where you run the code from.

## Requirements

To run the example code, it is assumed that the following R packages are installed on your system: 
  - R.Matlab
  - ggplot2
  - matrixcalc
