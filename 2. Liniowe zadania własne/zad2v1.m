clear;clc;
A = [4 2 0 0; 1 4 1 0; 0 1 4 1; 0 0 2 4];
n = length(A);
v = ones(n);

for k=1:10
    w=A*v;
    lam=(v'*w)/(v'*v);
    [y,i]=max(abs(w));
    alpha=w(i);
    v=w/alpha;
end
