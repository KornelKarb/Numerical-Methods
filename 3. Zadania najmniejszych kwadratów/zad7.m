clear;clc;

N1=10;
N2=-10;
N = [1;10;20;30;40;50;80;100];
l = [];
l2 =[];
for i = 1:8
c = round(linspace(0, N(i), N(i)));
r = round(linspace(0, -1*N(i), N(i)));
A = toeplitz(c,r);
x1 = round(linspace(0, N(i), N(i)));
x2 = rand(1, N(i));
b1 = A*x1';
b2 = A*x2';

alpha = 1e-4;
eps = 1e-6;
k = 10;

l(i) = norm(tikhonov(A,b1,alpha,eps));
l2(i) = norm(tikhonov(A,b2,k, eps));
end
condn1 = cond(A);
rankn1 = rank(A);
% 
tsvd1 = tsvd(A,b1,k, eps);
tikhonov1 = tikhonov(A,b1,alpha,eps);
% 
% SollerTsvdA = norm(tsvd1 -x1');
% ResidualTsvdA = norm(b1-A.*tsvd1 );
% SollerTikhonovA = norm(tikhonov1-x1');
% ResidualTikhonovA = norm(b1-A*tikhonov1);
% 
% tsvd2 = tsvd(A,b2,k, eps);
% tikhonov2 = tikhonov(A,b2,alpha,eps);
% 
% SollerTsvdB = norm(tsvd2 -x2');
% ResidualTsvdB = norm(b2-A.*tsvd2 );
% SollerTikhonovB = norm(tikhonov2-x2');
% ResidualTikhonovB = norm(b2-A*tikhonov2);

% y = [1e-7 ; 1e-5; 1e-4; 1e-3; 1e-2; 1e-1; 1; 10];
% l = [];
% l2 =[];
% l3 =[];
% l4=[];
% y2 = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
% for j = 1:8
%     l3(j) = norm(tikhonov(A, b2, y(j), eps));
%     l4(j) = norm(b2-A.*(tikhonov(A, b2, y(j), eps)));
% end
% for j = 1:k
% l(j) = norm(tsvd(A,b2,y2(j), eps));
% l2(j) = norm(b2-A.*(tsvd(A,b2,y2(j), eps)));
% end
figure
% semilogx(l2,l,'-o','LineWidth',1);
% loglog(l3,l4,'-o','LineWidth',1);
% loglog(l2,l,'-o','LineWidth',1);
plot(N,l2,'-o','LineWidth',1);
grid on
% plot(y,l);
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

