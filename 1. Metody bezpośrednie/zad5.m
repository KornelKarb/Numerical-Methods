clear;clc;
A = [1 2 3 4; -1 1 2 1; 0 2 1 3; 0 0 1 1];
[n,m] = size(A);
b = ones(n);
[L,U] = lu(A);
x = solve(A, b, L, U);

function [L,U] = lu(A)
m = size(A);
L = eye(m);
U = A;

for j=1:m-1
    for i=j+1:m
        L(i,j)=U(i,j)/U(j,j);
        U(i,j:m)=U(i,j:m)- L(i,j)*U(j,j:m);
    end
end
end

function x = solve(A, b, L)
n = size(b);
y = zeros(n);
x = zeros(n);

x(1)=b(1)/L(1,1);
for i=2:n
     x(i) = (b(i)-A(i,1:i-1)*x(1:i-1)')/A(i,i);
end
end