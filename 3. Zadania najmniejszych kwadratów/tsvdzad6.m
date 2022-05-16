clear;clc;
A = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
xt = [1 2 3]';
b = A*xt;
eps = 1e-6;
k = 3;
alpha = 1e-4;
x1 = tsvd(A,b,3, eps);
x12 = tsvd(A,b,2, eps);



function x = tsvd(A, b, k, eps)
[U, S, V] = svd(A);
n = size(S);
x = zeros(n);
for j = 1:k
    if S(j) >= eps
        x = x + (dot(U(:,j),b)\S(j))*V(j)';        
    end
end
end