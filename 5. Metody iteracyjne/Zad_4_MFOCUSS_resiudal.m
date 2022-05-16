clear; clc;
A = [2 3 -1 10 21 44 -9 1 -1; 1 2 2 8 15 35 8 -3 1; 3 1 1 6 16 53 -7 2 2];
b = [118 77 129]';

alpha = 1e-4;
p = 0;
lambda = 10e-3;
eps = 10e-3;
[m,n] = size(A);
X = rand(n,1);
p_kowergencji = 1;
k = 0;
X_poprz = X;
    while p_kowergencji > eps
        X_poprz = X;      
        w = sqrt(sum(abs(X).^2,2));
        W = diag(w.^(1-p/2));
        A = A*W;
        X = W*A'*inv(A*A'+lambda*eye(m))*b;
        p_kowergencji = sum(abs(X-X_poprz));
        k = k+1;
    end

x2 = tikhonov(A, b, alpha, eps);

SollerTikhonov = norm(x2-X');
ResidualTikhonov = norm(b-A*x2);


function x_przybl = tikhonov(A, b, alpha, eps)
m = length(A(1,:));
L = eye(m);
x_przybl = zeros(m,1);
error = 1;
    while error > eps
        x_poprz = x_przybl;
        x_przybl = x_przybl + (A'*A+alpha.^2*L)\A'*(b-A*x_przybl);
        error = norm(x_przybl-x_poprz);
    end
end

