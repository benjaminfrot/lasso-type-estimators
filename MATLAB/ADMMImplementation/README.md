An implementation of https://arxiv.org/abs/1512.06412 . We fit the model using the Alternative Direction Method of Multipliers (ADMM).

## Running the example code

From the MATLAB console, type:
```
run_example
```

The example script loads the toy dataset and computes solutions for a few values of lambda.

The function fit_many_lscggm.m is designed to be used as a standalone script, for example on a HPC cluster. It takes a single argument: the path to a .mat file which contains the data and the value of the tuning parameters to be used. It computes the solutions for the given list of tuning parameters and then outputs the results to a .mat file. See run_example.m for an example.

It is also possible to directly call fit_lscggm_with_split_bregman directly from MATLAB. Please see run_example.m and fit_many_lscggm.m for how to use these functions. In particular, the number of iterations, tolerance etc... can all be passed as arguments to fit_lscggm_split_bregman .

It is common that HPC clusters do not have MATLAB installed on every single node. Users are encouraged to use deploytool in order to create a standalone executable which can be executed even if MATLAB is not available (https://www.mathworks.com/help/compiler/deploytool.html). 

Questions? : benjamin dot frot at stat dot math dot ethz dot ch .
