clear;clc;
A = [4 2 0 0; 1 4 1 0; 0 1 4 1; 0 0 2 4];
numitr = 1000;
ep = 1.5;
[m,n]=size(A);
wynik = eigs(A);
y=ones(n,1);        
for k = 1 :  numitr
    iter = k; 
    v = y/norm(y,2);
    y = A\v;    
    th =v'*y; 
    if norm(y-th.*v,2) < ep*abs(th)
        break;
    end
end
  x = y/th;
