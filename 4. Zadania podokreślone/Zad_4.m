clear; clc;
A1 = [1 2 2 3 1; 2 4 4 6 2; 3 6 6 9 6;1 2 4 5 3];
A = [1 2 2 3 1; 0 4 4 6 2; 3 6 6 9 6;1 2 4 5 3];
X1 = [1 0 1 1 0]';
b = A*X1;

p = 1;
lambda = 10e-5;
X = mfocuss(A,b,lambda,p);

normalny = norm(X1 - X);
residualny = norm(b - A*X);

function X = mfocuss(A, b, lambda, p)

eps = 10e-5;
[m,n] = size(A);
X = rand(n,1);
p_kowergencji = 1;
k = 0;
X_poprz = X;
    while p_kowergencji > eps
        w = sqrt(sum(abs(X).^2,2));
        W = diag(w.^(1-p/2));
        A = A*W;
        %X = W*A'*(inv(A*A'+lambda*eye(m)).*b);
        X = W*A'*inv(A*A'+lambda*eye(m))*b;
        X_poprz = X;
        p_kowergencji = sum(abs(X-X_poprz));
        k = k+1;
    end

end