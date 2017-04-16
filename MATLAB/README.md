Implementations of the LSCGGM estimator described in (https://arxiv.org/abs/1512.06412).

This directory contains:
  - An implementation relying on the Alternative Direction Method Multipliers (in ADMMImplementation).
  - An implementation relying on Semi-Definite programming (in SDPImplementation).
  - R Code to generate data from the graphical models described in Section 5 of https://arxiv.org/abs/1512.06412 .

## Generating simulated data

You first need to install the R packages R.matlab and MASS, which are both available on CRAN:
```
install.packages("MASS")
install.packages("R.matlab")
```

The parameters of the R Script are: sample_size log2(rank of LZX) log2(rank of LX) outputfilename .
For example, running the following command will generate a dataset with 4 latent variables, in which inputs and outputs are in a one-to-one correspondence and save it to test_data.mat.
```
R --no-save --args 3000 5 2 ../test_data.mat < ./generate_example_data.R 
```
