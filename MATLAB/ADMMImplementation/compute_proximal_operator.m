function [X] = compute_proximal_operator(X, p, q, rho, tol, maxiter)

P = X * 0;
Q = X * 0;

for i=1:maxiter
    Y = prox_g(X + P, rho);
    P = X + P - Y;
    oX = X;
    X = prox_f(Y + Q, p, q);
    Q = Y + Q - X;
    if norm(oX, 2) < 1e-06
        diff = norm(X - oX, 2);
    else
        diff = norm(X - oX, 2) / norm(oX, 2);
    end
    if diff < tol
        break;
    end
end
end

function[A] = prox_g(X, rho)
    [U, S, V] = svd(X);
    S = S - rho;
    S(S < 0) = 0;
    A = U * S * V';
end

function [A] = prox_f(X, p, q)
    A = X(1:p,:);
    [V D] = eig(0.5 * (A + A'));
    D(D < 0) = 0;
    A = [V * D * V'; X((p+1):end,:)];
end