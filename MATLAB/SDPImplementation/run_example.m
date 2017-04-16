clear;

%% Uncomment this line to generate toy data. We use R (this requires the packages R.matlab and MASS)
% Make sure you run this script from the directory containing this script.
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

%% Fit the model using SDPT3 and the SDP implementation
% TO RUN THIS YOU NEED both YALMIP and SDPT3 to be in your path
% yalmip is available here: https://github.com/yalmip/YALMIP.git
% SDPT3 : https://github.com/sqlp/sdpt3.git
% Here we assume that both folders have been downloaded in the current
% directory
addpath(genpath('./YALMIP'))
cd sdpt3
install_sdpt3;
cd ..

% The function fit_many_lscggm fits the model for one or multiple values
% of the tuning parameters (lambda1, lambda2). It is designed to be easily used
% as part of a standalone executable deployed on a cluster with deploytool.
% For that reason, we write everything to a .mat file which is given as
% argument to the function. Of course, it is also possible to call
% fit_lscggm_with_SDPT3.m directly.
% Note that SDPT3 is *slow* and should not be used for large scale problems
% Use logdetPPA instead.
% Here, we compute a path.
gamma = 0.2;
lambda = [0.01 0.05 0.08];
% One pair (Lambda1, Lambda2) per line.
% Lambda1 is the penalty on ||S||_1, Lambda2 is the penalty on ||L||_\ast.
tuning_parameters = [lambda * gamma; lambda * (1 - gamma)]';

% Where we want the results to be written
output_filename = 'output.mat';
% Save the input data to input.mat
save('input.mat', 'output_filename', 'tuning_parameters', 'X', 'Z');
% Fit the model
tic
fit_many_lscggm('input.mat');
toc
% Load the results
load('output.mat')

%% Visualise the results.

%% Look at the sparsity patterns
% The top half is the estimate of SX. The bottom half is the estimate of
% SZX.
% With this solver, the estimate is not truly sparse: we threshold
spy(squeeze(abs(Ss(1,:,:))>1e-06)); % For first value of lambda
spy(squeeze(abs(Ss(2,:,:))>1e-06)); % For first value of lambda % For the the second value
spy(squeeze(abs(Ss(3,:,:))>1e-06)); % For first value of lambda % Third value

%% Look at the low-rank matrix
LX = squeeze(Ls(1,1:32,:)); % LX for the third value of Lambda
LZX = squeeze(Ls(1,33:end,:)); % LZX. The lower part of the matrix
%%
imagesc(LX); % Visualise it
%%
plot(svd(LX)) % Look at its spectrum