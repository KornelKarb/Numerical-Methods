clear;clc;
A = [-2,1,2,1; 2,-1,2,1; 2,3,-4,5; 2,3,0,-1];
% clear;clc;
% tic
% C = rand(200);
% A = kron(eye(30),C'*C);
% size(A);
% x= randn(6000,1);
% b=A*x;
% xr=A\b;
[Q,R]=gs(A);
%[W,T] = qr(A);
% toc
function [Q,R] = gs(A)
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(m,n);
    for j = 1:n
        Q(:,j) = A(:,j);
        for i = 1:j-1
            R(i,j) = Q(:,i)'*Q(:,j); %Tablica Q*Qt
            Q(:,j) = Q(:,j) - R(i,j)*Q(:,i);
        end
        R(j,j) = norm(Q(:,j))';
        Q(:,j) = Q(:,j)/R(j,j);
    end
end

