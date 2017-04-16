function [params, diffs] = fit_lscggm_with_split_bregman(Z, X, lambda1, lambda2, options, params)

maxiter = options.maxiter;
tol = options.tol;
nesterov_tol = options.nesterov_tol;
prox_tol = options.prox_tol;
prox_maxiter = options.prox_maxiter;

mu = 0.01;
p = size(X, 2);
q = size(Z, 2);
S = params.S;
L = params.L;
SX = S(1:p,:);
SZX = S((1+p):end, :);
LX = L(1:p, :);
LZX = L((1+p):end,:);
Lambda = params.Lambda;
diffs = [];
Theta	= struct;
Theta.yy = eye(p);
Theta.xy = zeros(q, p);

for i=1:maxiter
    % The first step : Update R
    [Theta, obj] = solve_step_1_nesterov(Z, X, SX, SZX, LX, LZX, Lambda, mu, nesterov_tol, Theta);
    RX = Theta.yy;
    RZX = Theta.xy;
    R = [RX; RZX];
    
    % Step 2: Proximal operator of the l1-norm: Prox(g, mu tau, S + G);
    A = (R + L + Lambda / mu);
    V = abs(A) - lambda1 / mu;
    V(V <= 0) = 0;
    S = sign(A) .* V;
    diff = norm(SX - S(1:p,:), 'fro') / norm(S,'fro');
    diffs = [diffs, diff];
    SX = S(1:p,:);
    SZX= S((p+1):end, :);
    
    % Step 3: Update L
    A = (S - R - Lambda / mu);
    L = compute_proximal_operator(A, p, q, lambda2 / mu, prox_tol, prox_maxiter);
    LX = L(1:p,:);
    LZX = L((p+1):end, :);
    
    % Step 4
    Lambda = Lambda + (R - S + L) * mu;
    
    params.S = S;
    params.L = L;
    params.R = R;
    params.Lambda = Lambda;
    
    if(diff < tol)
        break;
    end;
end;
end

