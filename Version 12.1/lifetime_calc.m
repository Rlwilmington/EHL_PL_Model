par = X;
A = par(1);
n = 1e13*par(2); %already converted to 1/m^2
T = par(3);
E_BGR = par(4);
alpha = par(5); %Lorentzian 'gamma' electrons -> alpha*(E-Ef)^2 + beta
beta = par(6);

%%
degen_e_K = 2; %band degeneracy
degen_h_K = 2;
degen_h_G = 1;

m0_oldunits = 1e6*0.510998950; %eV/c^2
c = 299792458; %m/s
m0 = (1/c^2)*m0_oldunits;
hbar = 6.582119569e-16; %eV*s
k_b = 8.617333262145e-5; %eV/K

E_opt = 1.89; %eV
E_bind = 0.421; %eV, note that exciton binding energy estimates vary 0.4-0.8eV (0.53 suggested by Robert Younts, 0.42 by Bataller et al)


const = [E_opt, E_bind, k_b, hbar, m0]; %constants necessary for fit model construction
%%

%alpha = 0.01;

Ef_e = F1(T,n,const,g_e_K);
Ef_h = F2(T,n,const,g_h_K,g_h_G,E_off_h);

gam_e = @(E) alpha*(E - Ef_e).^2 + beta;

gam_h = @(E) alpha*(E - Ef_h).^2 + beta;

N = 1000;

y1 = zeros(N,1);
y2 = zeros(N,1);

x = linspace(-0.1,0.1,N);

for i=1:N
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