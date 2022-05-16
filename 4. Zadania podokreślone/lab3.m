% clear; clc;
A = [1 2 2 3 1; 0 4 4 6 2; 3 6 6 9 6; 1 2 4 5 3];
x = max(0, rand(5,100));
c = length(x(1,:));
t = linspace(1,c,c);
bs = A*x;
figure
subplot(5,1,1);
plot(t,x(1,:));

subplot(5,1,2);
plot(t,x(2,:));

subplot(5,1,3);
plot(t,x(3,:));

subplot(5,1,4);
plot(t,x(4,:));

subplot(5,1,5);
plot(t,x(4,:));

o = linspace(1,5,5);

% q = zeros(5,5);
% z = zeros(5,5);
% u = zeros(5,5);
% resuf = zeros(5,1);
% soluf = zeros(5,1);
% resum = zeros(5,1);
% solum = zeros(5,1);

% for k =1:5
% xs = x(:,50+k);
% b = bs(:,50+k);
% foc= focuss(A, b, 0.001,0.1);
% mfoc = mfocuss(A, b, 1e-7,0.1);
% q(k,:) = foc';
% z(k,:) = mfoc';
% u(k,:) = xs';
% soluf(k) = norm(foc - xs);
% solum(k) = norm(mfoc - xs);
% resuf(k) = norm(b - A*foc);
% resum(k) = norm(b - A*mfoc);
% end
figure
subplot(2,1,1);
plot(o,q(1,:),o,u(1,:));
subplot(2,1,2);
plot(o,z(1,:),o,u(1,:));
figure
subplot(2,1,1);
plot(o,q(2,:),o,u(2,:));
subplot(2,1,2);
plot(o,z(2,:),o,u(2,:));
figure
subplot(2,1,1);
plot(o,q(3,:),o,u(3,:));
subplot(2,1,2);
plot(o,z(3,:),o,u(3,:));
figure
subplot(2,1,1);
plot(o,q(4,:),o,u(4,:));
subplot(2,1,2);
plot(o,z(4,:),o,u(4,:));
figure
subplot(2,1,1);
plot(o,q(5,:),o,u(5,:));
subplot(2,1,2);
plot(o,z(5,:),o,u(5,:));

function X = mfocuss(A, b, lambda, p)
eps = 10e-5;
[m,n] = size(A);
X = rand(n,1);
p_kowergencji = 1;
k = 0;
    while p_kowergencji > eps
        w = sqrt(sum(abs(X).^2,2));
        W = diag(w.^(1-p/2));
        A = A*W;
        %X = W*A'*(inv(A*A'+lambda*eye(m)).*b);
        X = W*A'*inv(A*A'+lambda*eye(m))*b;
        X_poprz = X;
        p_kowergencji = sum(abs(X-X_poprz));
        k = k+1;
    end

end

function X = focuss(A, b, lambda, p)
eps = 10e-5;
[m,n] = size(A);
%X = rand(n,1);
X = [1;1;1;1;1];
p_kowergencji = 1;
k = 0;
X_poprz = X;
    while p_kowergencji > eps
        W = diag(abs(X).^(1-p/2));
        X = W^2* A'*((A*W^2*A'+lambda*eye(m))\b);
        p_kowergencji = sum(abs(X-X_poprz));
        X_poprz = X;
        k = k+1;
    end
end