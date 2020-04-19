par = X;
A = par(1);
n = 1e13*par(2); %already converted to 1/m^2
T = par(3);
E_BGR = par(4);
alpha = par(5); %Lorentzian 'gamma' electrons -> alpha*(E-Ef)^2 + beta
beta = par(6);

%alpha = 0.01;

Ef_e = F1(T,n);
Ef_h = F2(T,n);

gam_e = @(E) alpha*(E - Ef_e).^2 + beta;

gam_h = @(E) alpha*(E - Ef_h).^2 + beta;

N = 1000;

y1 = zeros(N,1);
y2 = zeros(N,1);

x = linspace(-0.1,0.1,N);

for i=1:N;
    E = x(i);
    y1(i) = hbar/gam_e(E);
    y2(i) = hbar/gam_h(E);
end

hold on
plot(x,y1)
plot(x,y2)
legend('Electrons','Holes')
%ylim([1e-14, 1.3e-14])
hold off