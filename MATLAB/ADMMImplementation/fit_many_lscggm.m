function [] = fit_many_lscggm( input_filename )
% Forces MATLAB to use a single thread so that it can used as a standalone executable on a cluster
%maxNumCompThreads(1);

% Default settings.
maxiter = 500;
tol = 1e-04;
nesterov_tol = 1e-10;
prox_tol = 1e-10;
prox_maxiter = 100;

% Load the data
load(input_filename);
% Setup the options
options.maxiter = maxiter;
options.tol = tol;
options.nesterov_tol = nesterov_tol;
options.prox_tol = prox_tol;
options.prox_maxiter = prox_maxiter;
p = size(X,2);
q = size(Z,2);

% Initialise the parameters
SX = rand(p,p);
SX = SX * SX';
LX = zeros(p, p);
LZX = zeros(q, p);
SZX = zeros(q, p);
S = [SX; SZX];
L = [LX; LZX];
RX = SX - LX;
RZX = SZX - LZX;
R = [RX; RZX];
Lambda = zeros(p + q,p);
init = struct;
init.S = S;
init.L = L;
init.R = R;
init.Lambda = Lambda;

% Fit the model for all the values of the tuning parameters
% using the output of the previous fit as a warm start.
n_t = size(tuning_parameters,1);
Ss = zeros(n_t, p + q, p);
Ls = zeros(n_t, p + q, p);
Rs = zeros(n_t, p + q, p);
Lambdas = zeros(n_t, p + q, p);
Diffs = cell(n_t);
for i=1:n_t
    l1 = tuning_parameters(i,1);
    l2 = tuning_parameters(i,2);
    [l1 l2]
    [init, diffs] = fit_lscggm_with_split_bregman(Z, X, l1, l2, options, init);
    Ss(i,:,:) = init.S;
    Ls(i,:,:) = init.L;
    Rs(i,:,:) = init.R;
    Lambdas(i,:,:) = init.Lambda;
    Diffs{i} = diffs;
    save(output_filename, 'Ss', 'Ls', 'Rs', 'Lambdas', 'Diffs', 'tuning_parameters');
end;
end
