clear;clc;
A = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
xt = [1 2 3]';
b = A*xt;
eps = 1e-6;
k = 3;
alpha = 1e-4;
x1 = tsvd(A,b,k, eps);
x2 = tikhonov(A,b,alpha,eps);

SollerTsvd = norm(x1-xt');
ResidualTsvd = norm(b-A.*x1);

SollerTikhonov = norm(x2-xt');
ResidualTikhonov = norm(b-A*x2);

y = [1e-7 ; 1e-5; 1e-4; 1e-3; 1e-2; 1e-1; 1; 10];
l = [];
l2 =[];
l3 =[];
l4=[];
y2 = [1; 2; 3;];% 4; 5; 6; 7; 8; 9; 10];
for j = 1:8
    l3(j) = norm(tikhonov(A, b, y(j), eps));
    l4(j) = norm(b-A*(tikhonov(A, b, y(j), eps)));
end
for j = 1:k
l(j) = norm(tsvd(A,b,y2(j), eps));
l2(j) = norm(b-A.*(tsvd(A,b,y2(j), eps)));
end
figure
% semilogx(l2,l,'-o','LineWidth',1);
loglog(l4,l3,'-o','LineWidth',1);
% plot(y2,l);
grid on
%plot(y,l);
% plot(y,l);

function x = tsvd(A, b, k, eps)
[U, S, V] = svd(A);
n = size(S);
x = zeros(n);
S1 = diag(S);
for j = 1:k
    if S1(j) >= eps
        x = x + (dot(U(:,j)',b)/(S1(j)))*V(:,j)';        
    end
end
end

function x_przybl = tikhonov(A, b, alpha, eps)
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




% function y = tikhonov(A, b, k, alpha)
% 
% [U, S, V] = svd(A);
% n = size(S);
% y = zeros(n);
% for j = 1:k 
%     y = y + dot(U(:,j), b)*S(j)/(S(j)^2+alpha)*V(:,j)';
% end
% end