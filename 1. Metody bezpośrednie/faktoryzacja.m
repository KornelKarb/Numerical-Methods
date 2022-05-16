clear;clc;
%A = [1 2 3 4; -1 1 2 1; 0 2 1 3; 0 0 1 1];
%[n,m] = size(A);
% b = ones(n);
tic

%[T,W] = lu(A);
C = rand(200);
A = kron(eye(30),C'*C);
size(A);
x= randn(6000,1);
b=A*x;
xr=A\b;
[L,U]= dl(A);
%[x] = LUsolve(A,b,L,U);
toc
function [L,U] = dl(A)
[n,m] = size(A);
U = zeros(n);
L = zeros(n);
% %LU faktoryzacja metoda Dootlittle%
for j=1:m
    L(j,j)=1;
end
for j=1:m
    U(1,j)=A(1,j);
end
for i=2:m
    for j=1:m
        for k=1:i-1
            s1=0;
            if k==1
                s1=0;
            else
            for p=1:k-1
                s1=s1+L(i,p)*U(p,k); %wyznacznie sumy iloczynów wyrazów z L i U
            end
            end
            L(i,k)=(A(i,k)-s1)/U(k,k); %Obliczenie L(i,k)
         end
         for k=i:m
             s2=0;
           for p=1:i-1
               s2=s2+L(i,p)*U(p,k); %wyznacznie sumy iloczynów wyrazów z L i U
           end
           U(i,k)=A(i,k)-s2; %Obliczenie U(i,k)
         end
    end
end
end
% 
function [x] = LUsolve(A, b, L, U)
[n,m] = size(b);
y = zeros(n);
x = zeros(n);
for i=1:n
    y(i) = (b(i) - L(i,i)'*y(i))/L(i,i);
end
for i=n:1
    x(i)=(y(i)-U(i,i+1)'*x(i+1))/U(i,i);
end
end