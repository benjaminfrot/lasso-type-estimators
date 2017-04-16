An implementation of Equation (4.2) of https://arxiv.org/abs/1512.06412 . This is a semi-definite programming formulation
of the original problem.

## Requirements
We rely on two software packages: YALMIP and SDPT3. They must first be downloaded to this directory:
```
git clone https://github.com/yalmip/YALMIP
git clone https://github.com/sqlp/sdpt3.git
```

Please note that SDPT3 is convenient for prototyping but very slow and memory consuming. Consider using LogdetPPA for larger problems.

## Running the example code

From the MATLAB console, type:
```
run_example
```

The example script loads the toy dataset and computes solutions for a few values of lambda.

The function fit_many_lscggm.m is designed to be used as a standalone script, for example on a HPC cluster. It takes a single argument: the path to a .mat file which contains the data and the value of the tuning parameters to be used. It computes the solutions for the given list of tuning parameters and then outputs the results to a .mat file. See run_example.m for an example.

It is common that HPC clusters do not have MATLAB installed on every single node. Users are encouraged to use deploytool in order to create a standalone executable which can be executed even if MATLAB is not available (https://www.mathworks.com/help/compiler/deploytool.html). 
