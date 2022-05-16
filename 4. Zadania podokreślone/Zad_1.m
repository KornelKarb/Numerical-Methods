clear; clc;
% A = [1, 3, 1; -1, -2, 1; 3, 7, -1];
% X = [-2; 1; 0];
A1 = [1 2 2 3 1; 2 4 4 6 2; 3 6 6 9 6;1 2 4 5 3];
A2 = [1 2 2 3 1; 0 4 4 6 2; 3 6 6 9 6;1 2 4 5 3];
X1 = [1 0 1 1 0]';
b1 = A1*X1;
b2 = A2*X1;
p = 0;
lambda = 10e-5;
% X = focuss(A, b, lambda, p);
% normalny = norm(X1 - X);
% residualny = norm(b - A*X);

pwr = [0.0001,0.1,0.2,0.3,0.5,0.7,0.8,1];

xs = zeros(length(pwr));
ys = zeros(length(pwr));
xs2 = zeros(length(pwr));
ys2 = zeros(length(pwr));

for i = 1:length(pwr)
    c = focuss(A1, b1,lambda, pwr(i));
    ys(i) = norm(b1-A1*c);
    xs(i) = pwr(i);
    c1 = focuss(A2, b2,lambda, pwr(i));
    ys2(i) = norm(b1-A1*c1);
    xs2(i) = pwr(i);
end
figure
%plot(xs,ys);
plot(xs2,ys2);
grid on


function X = focuss(A, b, lambda, p)
eps = 10e-5;
[m,n] = size(A);
%X = rand(n,1);
X = [1;1;1;1;1];
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
end

