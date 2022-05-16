clear; clc;

% % % % % % Zadanie 1
% A = [1,2;4,2; -1,1];
% b = [4;12;1];
% c = [-1,-1];
%x =linspace(0,5,100);
% y1 = -0.5*x +2;
% y2 = -2*x+6;
% y3 = x+1;
% figure(1);clf;
% plot(x, y1);
% hold on;
% plot(x,y2);
% plot(x, y3);
% grid on;

% % % % % Zadanie 2
% A = [0.2,0.1;0.2,0.3];
% b = [60;120];
% c = [-0.03,-0.02];
% x =linspace(0,500,100);
% y1 = -0.5*x + (60/0.2);
% y2 = -1*(2/3)*x +(120/0.3);
% figure(1);clf;
% plot(x, y1);
% hold on;
% plot(y2,x);
% grid on;

% % % % % Zadanie 3
A = [0.5,2,2;1,1,1];
b = [7000;5000];
c = [-0.08,-0.01,-0.25];
x =linspace(0,10000,5000);
% y1 = (-2/0.5)*x +(-2/0.5)*x + (7000/0.5);
% y2 = -x -x +5000;
% y3 = (-0.08/0.25)*x + (-0.1/0.25)*x +(100/0.25);
% figure(1);clf;
% plot3(y1, y2,y3);
% hold on;
% grid on;

% % % % % Zadanie 4
% A = [1,0;0,1;1,-2];
% b = [6.4*10^6;3*10^6;0];
% c = [-1.9,-1.5];


% % % % % Zadanie 5
% A = [1,1,1;1,0,0;0,1,-3];
% b = [12000;2000;0];
% c = [-1.12,-1.07,-1.08];


% % % % % Zadanie 6
% A = [2,1;1,3;1,0];
% b = [40;45;12];
% c = [-30,-25];
% x =linspace(0,15,2000);
% y1 = 40 - 2*x;
% y2 = 45/3 - (1/3)*x;

[w,T] = simplex(A,b,c);
n = linprog(c,A,b);

% figure(1);clf;
% plot(x, y1);
% hold on;
% plot(x,y2);
% plot(x, y3);
% grid on;

function [x,T] = simplex(A,b,c)
[m,n] = size(A);
diag = eye(m);
T = cat(2, A, diag);
T = [T b];
ct = cat(2,c,zeros(1,m+1));
T = vertcat(T,ct);
disp(T);

[row,col] = size(T);


    for i = 1:col
     val = min(T(row,:));
     pivot_row = 0;
        for z=1:col
            if T(row,z) == val
            col_inx = z;
        break;
            end
        end
        
         if val < 0
             pivot_row = findPivot(T, col_inx);
         else
             x = zeros(n);
             for it = 1:n
                 piv = max(T(:,it));
                 spt = find(T(:,it) == piv);
                 if T(row, it) == 0
                     x(it) = T(spt,col);
                 else
                     x(it) = 0;
                 end
             end
         end
        
         if pivot_row > 0
            T(pivot_row,:) = T(pivot_row,:)/T(pivot_row, col_inx);
            for k = 1:m+1
                if k ~= pivot_row
                    %T(k,:) = T(k,:)-T(k,col_inx)*T(pivot_row,:);
                    T(k,:) = T(k,:)-sign(T(k, col_inx))*abs(T(k, col_inx))*T(pivot_row,:);
                end
            end
         end
         disp(T);
    end
    
end

function pivot = findPivot(T, column)
[m,n] = size(T);
pivot = 0;
min = inf;
for k = 1:m
    if T(k,column) > 0
        tmp = T(k, end)/T(k,column);
        if tmp < min
            min = tmp;
            pivot = k;
        end
    end
end
end

