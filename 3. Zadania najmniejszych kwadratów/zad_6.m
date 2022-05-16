clear;clc;

A = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
x = [1; 2; 3];
b = A*x;
alpha = 1e-6;
eps = 1e-6;

condn1 = cond(A);
[U, S, V] = svd(A);
condn2 = S(1)/S(3,3);

test_przybl = tikhonova(A, b, alpha, eps);

SollerTikhonov = norm(test_przybl-x');
ResidualTikhonov = norm(b'-A*test_przybl);

y = [1e-7 ; 1e-5; 1e-4; 1e-3; 1e-2; 1e-1; 1; 10];
l = [];
l2 = [];
y2 = [1; 2; 3];
for j = 1:8
    l(j) = norm(tikhonova(A, b, y(j), eps));
end

figure
semilogx(y,l);
grid on
% plot(y,l);



function x_przybl = tikhonova(A, b, alpha, eps)
m = length(A(1,:));
L = eye(m);
x_przybl = zeros (m,1);
error = 1;
    while error > eps
        x_poprz = x_przybl;
        x_przybl = x_przybl + (A'*A+alpha.^2*L)\A'*(b-A*x_przybl);
        error = norm(x_przybl-x_poprz);
    end
 
end
