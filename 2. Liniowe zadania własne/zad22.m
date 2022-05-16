clear;clc;
A = [4 2 0 0; 1 4 1 0; 0 1 4 1; 0 0 2 4];
[m,n]=size(A);
v = ones(n);
s=1.5;
As=A-s*eye(n);
for i=1:100
    w =As\v;
    lam=(v'*v)/(v'*w)+s;
    [y,i] = max(abs(w));
    alpha=w(i);
    v=w/alpha;
end