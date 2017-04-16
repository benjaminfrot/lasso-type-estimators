function [] = fit_many_lscggm( input_filename )
% Forces MATLAB to use a single thread so that it can used as a standalone executable on a cluster
maxNumCompThreads(1);

% Load the data
load(input_filename);
p = size(X,2);
q = size(Z,2);
Sigma = cov([Z X]);
%% Set up the problem with yalmip
addpath(genpath('./yalmip'));
addpath(genpath('./sdpt3'));

% Set up the problem
SX = sdpvar(p,p);
SZX = sdpvar(q,p);
LX = sdpvar(p,p);
LZX = sdpvar(q,p);
W = sdpvar(q,q);
S = [W (SZX - LZX); (SZX - LZX)' (SX - LX)];
L = [LZX; LX];
% Now introduce additional variables for the nuclear norm
W1 = sdpvar(p + q,p + q);
W2 = sdpvar(p,p);
C = [W1 L; L' W2];

Constraints = [S >= 0, SX - LX >= 0, LX >= 0, C >= 0];


%% Fit the model for all the values of the tuning parameters
% using the output of the previous fit as a warm start.
n_t = size(tuning_parameters,1);
Ss = zeros(n_t, p + q, p);
Ls = zeros(n_t, p + q, p);
for i=1:n_t
    l1 = tuning_parameters(i,1);
    l2 = tuning_parameters(i,2);
    [l1 l2]
    
    Objective = trace(S * Sigma) - logdet(SX - LX) + l1 * sum(sum(abs([SZX; SX]))) + l2 * 0.5 * (trace(W1) + trace(W2));
    options = sdpsettings('verbose',1,'solver','sdpt3', 'dualize', 0);
    sol = optimize(Constraints,Objective,options);
    SXv = value(SX);
    LXv = value(LX);
    SZXv = value(SZX);
    LZXv = value(LZX);
    Ss(i,:,:) = [SXv; SZXv];
    Ls(i,:,:) = [LXv; LZXv];
    save(output_filename, 'Ss', 'Ls', 'tuning_parameters');
end;
end