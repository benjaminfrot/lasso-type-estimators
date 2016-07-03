function [] = fit_one_with_SDPT3(input_filename)
%%
% Takes a path to a .mat file as input
% The .mat file must contain the following variables
% - use_covariance : logical. Indicates whether the covariance matrix
% Sigma is provided, or whether the raw data Z, X is given
% - if use_covariance=1, Sigma : a (p + q) x (p + q) positive
% semi-definite matrix. Columns/Rows 1 to q correspond to Z. The rest
% to X
% - if use_covariance=1, q : an integer. The number of inputs |Z|. (The
% number of outputs |X| is deduced from q and the size of Sigma)
% - if use_covariance=0, Z and X: n x q and n x p matrices containing
% the raw data. One sample per line. They will be standardised.
% - output_filename : name of file where to write the results
% - tuning_parameters : a k x 2 matrix. The list of tuning parameters
% to use. First column is the l1-penalty. Second column is the
% l2-penalty


%%
% maxNumCompThreads(1); % Used on clusters to prevent jobs from killed
% because they use more than one thread

%%
% Assumes that yalmip and sdpt3 are in this folder
addpath(genpath('./yalmip'));
addpath(genpath('./sdpt3'));
%%
% Load everything
load(input_filename);

%%
if use_covariance == 0
    q = size(Z,2);
    p = size(X,2);
    for i=1:p
        v = X(:,i);
        v = (v - mean(v)) / sqrt(var(v));
        X(:,i) = v;
    end;
    for i=1:q
        v = Z(:,i);
        v = (v - mean(v)) / sqrt(var(v));
        Z(:,i) = v;
    end;
    Sigma = cov([Z X]);
else
    p = size(Sigma,2) - q;
end;

%% The model
SX = sdpvar(p,p);
SZX = sdpvar(q,p);
LX = sdpvar(p,p);
LZX = sdpvar(q,p);
W = sdpvar(q,q);
Zc = sdpvar(p+q,p);
S = [W (SZX - LZX); (SZX - LZX)' (SX - LX)];
L = [LZX; LX];

% Now introduce additional variables for the nuclear norm
W1 = sdpvar(p + q,p + q);
W2 = sdpvar(p,p);
C = [W1 L; L' W2];

Constraints = [S >= 0, SX - LX >= 0, LX >= 0, C >= 0, [SZX; SX] <= Zc, [SZX; SX] >= -Zc];

n_t = size(tuning_parameters,1);
results = cell(n_t);
for i=1:n_t
    l1 = tuning_parameters(i,1);
    l2 = tuning_parameters(i,2);
    
    Objective = trace(S * Sigma) - logdet(SX - LX) + ...
        l1 * ones(1,p+q) * Zc * ones(p,1) + ...
        l2 * 0.5 * (trace(W1) + trace(W2));
        options = sdpsettings('verbose',1,'solver','sdpt3', 'dualize', 0, ...
            'sdpt3.gaptol', 1e-08,  'sdpt3.inftol', 1e-08);
    sol = optimize(Constraints,Objective,options);
    
    results{i} = struct;
    results{i}.sol = sol;
    results{i}.SX = value(SX);
    results{i}.SZX = value(SZX);
    results{i}.LX = value(LX);
    results{i}.LZX = value(LZX);
end;

save(output_filename, 'results', 'tuning_parameters');
end