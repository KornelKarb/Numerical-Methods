clear;clc;
% A = [2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
% b = [ 0; 0; 0; 5];
% 
% A = [A b];
% [m,n] = size(A);
tic
C = rand(200);
A = kron(eye(30),C'*C);
[m,n]=size(A);
x= randn(6000,1);
b=A*x;
xr=A\b;
for k = 1:m
    MaxElem = 0;
    MaxWierszIndex = 1;
    for i=k:m
        IndexNow = i;
        ElemNow = abs(A(i,k));
        if ElemNow >= MaxElem
            MaxWierszIndex = IndexNow;
            MaxElem = ElemNow;
        end
    end
tempRow = A(k,:);
A(k,:) = A(MaxWierszIndex,:);
A(MaxWierszIndex,:) = tempRow;

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
end
