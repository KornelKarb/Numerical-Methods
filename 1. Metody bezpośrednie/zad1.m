clear;clc;
A = [2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
b = [ 0; 0; 0; 5];

A = [A b];
[m,n] = size(A);

%Eliminacja gaussa%
%Pętla eleminacji%
for i=1:m-1
    for j=i+1:m
        %A(i,:)=A(i,:)-A(j,:)*(A(i,j)/A(j,j));%
        Multipliers = A(j,i)/A(i,i);
        for l = 1:m+1
            A(j,l) = A(j,l) - Multipliers * A(i,l);
        end
    end
end
x=zeros(1,m);
%Wyznaczenie wartości dla zmiennych%
for i=m:-1:1
    x(i)=A(i,m+1);
    for j=i+1:m
        x(i)= x(i)-A(i,j)*x(j);
    end
    x(i)=x(i)/A(i,i);
end