clear; clc;
tic
% A = [1 1 1; 1 1 2; 1 2 2];
% b = [1; 2; 1];
% x_nom = A\b;

% % Zadanie 5
% n = 100;
% D1 = diag(diag(2*ones(n)));
% D2 = diag(diag(-1*ones(n-1)),1);
% D3 = diag(diag(-1*ones(n-1)),-1);
% A = D1+D2+D3;
% b = ones(n);
% x_nom = A\b;

% %zadanie6
C = rand(200);
A = kron(eye(10),C'*C);
size(A);
x= randn(2000,1);
b=A*x;
x_nom = A\b;

x0 = 0; %wartość początkowa
kmax = 100;
omega = 1.2;
alpha = 0.1111;
tsolve = linsolve(A,b);
% [kowerg, blad_resyd, zbiez] = gauss(A, b, x_nom, kmax, x0);
% [kowergj, blad_resydj, zbiezj] = jacobi(A, b, x_nom, kmax, x0);
% [x, res_err, zbiezsd] = sd(A, b, kmax, tsolve);
% [kowergk, blad_resydk, zbiezk] = kaczmarz(A, b, tsolve, kmax, omega, x0);
% [kowergs, blad_resyds, zbiezs] = SOR(A, b, tsolve, kmax, omega, x0);
[kowergl, blad_resydl, zbiezl] = landweber(A, b, tsolve, kmax, alpha, x0);

% figure(1); clf;
% title('Kowergencja');
% hold on;
% semilogy(kowergj);
% grid on;
% legend('Gauss');
% 
% figure(2); clf;
% title('Kowergencja');
% hold on;
% semilogy(zbiezj);
% grid on;
% legend('Gauss');
% 
% figure(3); clf;
% title('Kowergencja');
% hold on;
% semilogy(blad_resydj);
% grid on;
% legend('Gauss');
toc
% figure(1); clf;
% semilogy(kowergk);
% title('Kowergencja');
% hold on;
% semilogy(kowergj);
% semilogy(kowergs);
% semilogy(kowergl);
% semilogy(kowerg);
% grid on;
% legend('Kaczmarz','Jacobi','SOR','Landweber','Gauss');
% 
% figure(2);clf;
% semilogy(zbiezk);
% title('Błąd rozwiązania');
% hold on;
% semilogy(zbiezj);
% semilogy(zbiezs);
% semilogy(zbiezl);
% semilogy(zbiez);
% grid on;
% legend('Kaczmarz','Jacobi','SOR','Landweber','Gauss');
% 
% figure(3);clf;
% semilogy(blad_resydk);
% title('Błąd Residualny');
% hold on;
% semilogy(blad_resydj);
% semilogy(blad_resyds);
% semilogy(blad_resydl);
% semilogy(blad_resyd);
% grid on;
% legend('Kaczmarz','Jacobi','SOR','Landweber','Gauss');


% metoda Gaussa
function [kowerg, blad_resyd, zbiez] = gauss(A, b, x_nom, kmax, x0)
S = tril(A);
T = -triu(A,1);
x = x0;
for k = 1:kmax
    G = S\T;
    c = S\b;
    xprev = x;
    x = G*x + c;

    kowerg(k) = norm(x(:,1) - x_nom, 2);
    blad_resyd(k) = norm(b-A*x, 2);
    zbiez(k) = norm(x-xprev, 2);
end
end

% metoda Jacobi
function [kowergj, blad_resydj, zbiezj] = jacobi(A, b, x_nom, kmax, x0)
% x = x(:,1);
xj = x0;
Sj = diag(diag(A));
Tj = Sj - A;
x = x0;
for k = 1:kmax
    Gj = Sj\Tj;
    cj = Sj\b;
    xprevj = xj;
    xj = Gj*xj + cj;
    kowergj(k) = norm(xj(:,1) - x_nom, 2);
    blad_resydj(k) = norm(b-A*xj, 2);
    zbiezj(k) = norm(xj - xprevj, 2);
end
end
%Kaczmarz
function [kowergk, blad_resydk, zbiezk] = kaczmarz(A, b, tsolve, kmax, omega, x0)
m = size(A);
x = zeros(m);

for i = 1:kmax
    for j = 1:m
        x = x + omega*((b(j)-A(j).*x)/(norm(A(j))^2))*A(j);
    end
    blad_resydk(i) = norm(b-A.*x, 2);
    zbiezk(i) = norm(x - tsolve, 2);
    kowergk = x;
end
end
%SOR
function [kowergs, blad_resyds, zbiezs] = SOR(A, b, tsolve, kmax, omega, x0)
m = size(A);
L = A;
D = eye(m);
D = A*D;

for i = 1:m
    for j = 1:m
    if j > i
        L(i,j)=0;
    end
    end
end
U = A - L;
L = L - D;
x = x0;

for i = 1:kmax
    x = inv(L + omega*D).*(omega*b-(omega*U+(omega-1)*D).*x);
    blad_resyds(i) = norm(b-A*x, 2);
    zbiezs(i) = norm(x - tsolve, 2);
    kowergs(i) = norm(x);
end
end
%Landweber
function [kowergl, blad_resydl, zbiezl] = landweber(A, b, tsolve, kmax, alpha, x0)
m = size(A);
x = x0;

for i = 1:kmax
    x = x + alpha*A'.*(-A.*x+b);
    blad_resydl(i) = norm(b-A.*x, 2);
    zbiezl(i) = norm(x - tsolve, 2);
    kowergl(i) = norm(x);
end
end
%%Algorytm SD
function [x, res_err, zbiez] = sd(A, b, kmax, tsolve)
m = size(A);
x = ones(m);
res_err = [];
zbiez = [];

for i = 1:kmax
    r = b - (A.*x);
    f = (r'.*r)/(r'.*A.*r);
    x = x + f*r;
    res_err(i) = norm(b-A*x, 2);
    zbiez(i) = norm(x - tsolve, 2);
end
end