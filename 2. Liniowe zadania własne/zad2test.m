clear;clc;
A = [4 2 0 0; 1 4 1 0; 0 1 4 1; 0 0 2 4];
n = length(A);

z = ones(n);
q = z/sqrt(sum(z'*z));

for k = 1:5
    z = A*q;
    q=z/sqrt(sum(z'*z));
end
lam = 0;
m = 0;
for i = 1:n
  lam = lam + z(i)*q(i);
  m = m + lam*z(i);
end


