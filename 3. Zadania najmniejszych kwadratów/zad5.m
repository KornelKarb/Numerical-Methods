x = linspace(0, pi);
F = pi.^2 - x.^2;
f = F.';
n = 10;

G = zeros(length(f),n);

for i= 0:n
    G(:,i+1) = cos(i*x);
end

c = lsqr(G,f);

y = G*c;

figure
plot(x,f,x,y)
grid on