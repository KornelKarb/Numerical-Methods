clear; clc;
A = [1 2 2 3 1; 2 4 4 6 2; 3 6 6 9 6;1 2 4 5 3];
b = [1 0 1 1 0]';
p = 1;
lambda = 10e-5;
eps = 10e-5;
[m,n] = size(A);
X = rand(n,1);
p_kowergencji = 1;
k = 0;
X_poprz = X;
    while p_kowergencji > eps
        W = diag(abs(X).^(1-p/2));
        X = W^2* A'*((A*W^2*A'+lambda*eye(m))\b);
        p_kowergencji = sum(abs(X-X_poprz));
        X_poprz = X;
        k = k+1;
    end