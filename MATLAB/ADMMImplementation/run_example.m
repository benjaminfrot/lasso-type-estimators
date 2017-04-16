clear;

%% Uncomment this part to generate toy data. We use R (this requires the packages R.matlab and MASS)
% Make sure you run this script from the directory containing it.
% By default |X| = |Z| = 32
% The parameters are : sample_size log2(rank of LZX) log2(rank of LX)
% output_filename
% In this intance: SZX is sparse (2**5 = 32) and LX has rank 4 (2**2 = 4)
% ! R --no-save --args 3000 5 2 ../test_data.mat < ../generate_example_data.R 

%% Load some test_data
load ../test_data;
% SX is the true sparse matrix:
spy(SX);
% LXmle is the true latent matrix LX computed from 1e5 samples. It's rank
% is 2**u
% Likewise KZXmle is KZX = SZX - LZX computed from 1e5 samples. LZXmle =
% 2**z
% Sig_XH, Sig_Z, Sig_Z_XH are covariance matrices computed from 1e5
% samples. They also include the hidden variables.
% Finally Z and X: data we use to learn the model. ObsData = [Z X];

%% Fit the model using the Split-Bregman method
% The function fit_many_lscggm fits the model for one or multiple values
% of the tuning parameters (lambda1, lambda2). It is designed to be easily used
% as part of a standalone executable deployed on a cluster with deploytool.
% For that reason, we write everything to a .mat file which is given as
% argument to the function. Of course, it is also possible to call
% fit_lscggm_with_split_bregman.m directly.

% Here, we compute a path.
gamma = 0.3;
lambda = [0.01, 0.05, 0.1];
% One pair (Lambda1, Lambda2) per line.
% Lambda1 is the penalty on ||S||_1, Lambda2 is the penalty on ||L||_\ast.
tuning_parameters = [lambda * gamma; lambda * (1 -gamma)]';

% Where we want the results to be written
output_filename = 'output.mat';
% Maximum number of iterations. See function for more options. By default
% tolerance is 10-5 . 
maxiter = 1000;
% Save the input data to input.mat
save('input.mat', 'output_filename', 'tuning_parameters', 'X', 'Z', 'maxiter');
% Fit the model
tic
fit_many_lscggm('input.mat');
toc
% Load the results
load('output.mat')

%% Visualise the results.

%% Look at the sparsity patterns
% The top half is the estimate of SX. The bottom half is the estimate of
% SZX
spy(squeeze(Ss(1,:,:))); % For first value of lambda
spy(squeeze(Ss(2,:,:))); % For the the second value
spy(squeeze(Ss(3,:,:))); % Third value

%% Look at the low-rank matrix
LX = squeeze(Ls(2,1:32,:)); % LX for the third value of Lambda
LZX = squeeze(Ls(2,33:end,:)); % LZX. The lower part of the matrix
%%
imagesc(LX); % Visualise it
%%
plot(svd(LX)) % Look at its spectrum

%% Look at the convergence speed of the algorithm by plotting the change in the norm of the parameters between consecutive iterations
plot(Diffs{3});

%% Look at the sparsity of S along the path
plot(sum(sum(Ss~=0,3),2));

%% Look at the number of iterations required to reach convergence with a tolerance of 1-04, as a function of lambda
[nrows, ncols] = cellfun(@size, Diffs);
plot(ncols(:,1));