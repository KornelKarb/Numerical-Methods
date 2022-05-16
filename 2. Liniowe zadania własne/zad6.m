clc;clear;
A = [1 1; 0 1; 1 0];
% A=[1 2 2];
% A = [ 2 2 2 2; 1.7 0.1 -1.7 -0.1;0.6 1.8 -0.6 -1.8];

[n,m] = size(A);
U = eye(n);
V = eye(m);
tol = max(abs(A(:)))*1.e-15;
Arem = inf;
while Arem > tol
    [Qu,Ru] = qr(A);
    U = U*Qu;
    [Qv,Rv] = qr(Ru');
    V = V*Qv;
    A = Rv';
    Arem = norm(tril(A,-1),inf);
end
U = U(:,1:m);
S = triu(A(1:m,:));
V = V.*sign(diag(S)).';

