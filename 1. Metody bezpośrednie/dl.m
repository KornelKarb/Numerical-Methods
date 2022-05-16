tic
C = rand(200);
A = kron(eye(30),C'*C);
size(A);
x= randn(6000,1);
b=A*x;
xr=A\b;
[L,U]= dla(A);
toc
function [L,U] = dla(A)
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
            sum1=0;
            if k==1
                sum1=0;
            else
            for p=1:k-1
                sum1=sum1+L(i,p)*U(p,k);
            end
            end
            L(i,k)=(A(i,k)-s1)/U(k,k);
         end
         for k=i:m
             sum2=0;
           for p=1:i-1
               sum2=sum2+L(i,p)*U(p,k);
           end
           U(i,k)=A(i,k)-sum2;
         end
    end
end
end
