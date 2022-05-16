clear;clc;
tic
C = rand(200);
A = kron(eye(30),C'*C);
size(A);
x= randn(6000,1);
b=A*x;
xr=A\b;
[Q,R] = qr(A);
toc