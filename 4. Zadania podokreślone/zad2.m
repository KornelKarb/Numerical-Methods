clear; clc;
A = [1 2 2 3 1; 2 4 4 6 2; 3 6 6 9 6;1 2 4 5 3];
b = [1 0 1 1 0]';
lam = 0.000001;
[m,n] = size(A);
x = rand(n,1);
I = eye(m);
lastx = 10*x;
power = 2;
for k = 1:500
    lastx = x;
    w = sqrt(sum(abs(x).^2,1));
    W = (w.^(1-(power/2)))*eye(n);
    A = A*W;
    Q = A'.*inv(A.*A'+lam*I).*b;
    x = W.*Q;
    xs = x;
end