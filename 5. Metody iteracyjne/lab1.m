clear; clc;
A = [1 1 1; 1 1 2; 1 2 2];
b = [1; 2; 1];
x_nom = A\b;

x0 = 0; %wartość początkowa
kmax = 30;
tsolve = linsolve(A,b);

m = size(A);
x = ones(m);
res_err = [];
zbiez = [];

for i = 1:kmax
    r = b - (A*x);
    f = (r'*r)/(r'*A*r);
    x = x + f*r;
    res_err(i) = norm(b-A*x, 2);
    zbiez(i) = norm(x - tsolve, 2);
end
