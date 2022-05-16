clear; clc;
A = [1 3 1; -1 -2 1; 3 7 -1];
b = setvar(0, 3, 1);
lam = 0.000001;
tolerance = 1e-10;
x_star = [0 3 1];
x1 = focuss(A, b, lam, tolerance)
function b = setvar(x, y, z)
    a = x+3*y+z;
    b = -x-2*y+z;
    c = 3*x+7*y-z;
    b = [a,b,c];
end

function x = focuss(A, b, lam, power, tolerance)
s = size(A);
I = eye(s);
x = rand(1,1);
x = x(:,0);
lastx = 10*x;

while norm(lastx-x)/norm(x) > tolerance
    lastx = x;
    W = diag(x)^(1-(power/2));
    W2 = W^2;
    x = W2.*A'.*inv(A.*W2.*A+lam*I).*b;
end
end


