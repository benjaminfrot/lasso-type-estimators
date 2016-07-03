% Load the dataset
clear;
load example_dataset.mat;

output_filename = 'output.mat';
use_covariance = 0;
lambda = 0.01;
gamma = 0.5;
tuning_parameters = [lambda * gamma, lambda * (1 - gamma)];

% Write everything to a file called input.mat
save('input.mat', 'output_filename', 'Z', 'X', 'use_covariance', 'tuning_parameters');

% Fit with the SDPT3 implementation of LSCGGM
fit_lscggm_with_SDPT3('input.mat')

load output.mat
% Look at SZX
imagesc(abs(results{1}.SX) > 1e-06)

%% Now, let's do the same, but this let's compute the covariance matrix ourselves
clear;
load example_dataset.mat;

output_filename = 'output.mat';
use_covariance = 1;
q = size(Z,2);
p = size(X,2);
for i=1:q
    v = Z(:,i);
    v = (v - mean(v)) / sqrt(var(v));
    Z(:,i) = v;
end;
for i=1:p
    v = X(:,i);
    v = (v - mean(v)) / sqrt(var(v));
    X(:,i) = v;
end;
Sigma = cov([Z X]);

lambdas = [0.1, 0.01, 0.001];
gamma = 0.5;
tuning_parameters = [];
for i=1:length(lambdas)
       tuning_parameters = [tuning_parameters; [lambdas(i) * gamma, lambdas(i) * (1 - gamma)]];
end;
tuning_parameters
% Write everything to a file called input.mat
save('input.mat', 'output_filename', 'Sigma', 'q', 'use_covariance', 'tuning_parameters');

% Fit with the SDPT3 implementation of LSCGGM
fit_lscggm_with_SDPT3('input.mat')

load output.mat
% Look at SZX;
imagesc(abs(results{3}.SX) > 1e-06)
