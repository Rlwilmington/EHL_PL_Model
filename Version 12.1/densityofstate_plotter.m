%% Set Up Model Fit

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
X = [4.8806  613.3210    0.6497    8.6656    1.2369    0.4784   -0.0000   -0.0036];
par = X;

[f1, f2, f3, f4, f5] = StrainvTemp_MAC();

n = 1e13*par(1); %already converted to 1/m^2
T = par(2);
E_BGR = par(3);
alpha = par(4); %Lorentzian 'gamma' -> alpha*(E-Ef)^2 + beta
beta = par(5);
A = 1e-38*par(6);
C = par(7);
D = par(8);

step = 0.005;
xfit_g = (-5:step:5).'; %end at 4eV, though this choice should not make a difference

if X(2) > 294
    strain = polyval(f1,X(2));
else
    strain = polyval(f1,294);
end


disp('Sample Strain (Lattice Expansion)')
disp(strain)
E_off_h = polyval(f3,strain); %eV

%% Electrons/Holes Fermi Levels
m_e_K = polyval(f5,strain)*m0; %effective masses (highest fluence)
m_h_K = polyval(f4,strain)*m0;
m_h_G = polyval(f2,strain)*m0;
g_h_K = (degen_h_K*m_h_K)/(pi*hbar^2);
g_h_G = (degen_h_G*m_h_G)/(pi*hbar^2);
g_e_K = (degen_e_K*m_e_K)/(pi*hbar^2); %density of states

E_fermi_e = F1(X(2),X(1).*1e13,const,g_e_K);
E_fermi_h = F2(X(2),X(1).*1e13,const,g_h_K,g_h_G,E_off_h);






%% Plot

u1 = @(E) D_e(const,E,n,T,alpha,C,f1,f5);
v1 = @(E) D_h(const,E,n,T,beta,D,f1,f2,f3,f4);

u2 = zeros(length(xfit_g),1);
v2 = zeros(length(xfit_g),1);


for i = 1:length(xfit_g)
    
    E = xfit_g(i);
    u2(i) = D_e_noL(const,E,n,T,f1,f5);
    v2(i) = D_h_noL(const,E,n,T,f1,f2,f3,f4);
end

u1_vec = u1(xfit_g);
v1_vec = v1(xfit_g);

% xfit_g_e1 = E_opt + xfit_g; %flipud(xfit_g);
% xfit_g_h1 = flipud(xfit_g) + E_off_h;
% 
% xfit_g_e2 = E_opt + flipud(xfit_g); %flipud(xfit_g);
% xfit_g_h2 = xfit_g + E_off_h;

xfit_g_e1 = E_opt + flipud(xfit_g);
xfit_g_e2 = E_opt + flipud(xfit_g);
xfit_g_h1 = -E_off_h + xfit_g;
xfit_g_h2 = -E_off_h + xfit_g;

u1_vec(u1_vec < 1e15) = -1;
v1_vec(v1_vec < 1e15) = -1;
u2(u2 < 1e15) = -1;
v2(v2 < 1e15) = -1;

scale = 1;

u2 = scale*u2;
v2 = scale*v2;

%E_gap = E_opt + E_bind - E_off_h - E_BGR - E_thermal;

E_gap = E_opt - E_off_h;

f = figure('Position',[100 100 900 400]);
hold on
plot(E_gap + 2*E_fermi_e + flipud(xfit_g),u1_vec,'k-','LineWidth',2); %convolved electron
plot(E_gap + 2*E_fermi_e + flipud(xfit_g),u2,'k--','LineWidth',2); %raw electron
plot([E_gap + 2*E_fermi_e + -E_fermi_e, E_gap + 2*E_fermi_e + -E_fermi_e],[0 6e18],'r-','LineWidth',2) %Ef electron
plot([0, 0],[0, 2e18],'b:','LineWidth',2) %gap
plot(-E_off_h + xfit_g,v1_vec,'k-','LineWidth',2);
plot(-E_off_h + xfit_g,v2,'k--','LineWidth',2);
plot([-E_off_h - E_fermi_h,-E_off_h - E_fermi_h],[0 6e18],'-r','LineWidth',2)
plot([E_gap, E_gap],[0, 2e18],'b:','LineWidth',2)
plot([0, E_gap],[2e18, 2e18],'b:','LineWidth',2)
xlim([-0.5 2.5])
ylim([1e10 6e18])
legend('Lorentzian Convolved Densities of States','Densities of States (no broadening)','Quasi-Fermi Level','Energy Gap','Location','North')
xlabel('Energy (eV, referenced at K-VB)')
ylabel('Density of States (arb.)')
set(gca,'FontSize',14)
hold off

saveas(gcf,'fig_DOS.png')


% figure();
% hold on
% plot(xfit_g_e1,u1_vec)
% plot(xfit_g_h1,v1_vec)
% plot(xfit_g_e2,u2,'k--')
% plot(xfit_g_h2,v2,'k--')
% plot([-E_off_h + E_opt+E_fermi_e, -E_off_h + E_opt+E_fermi_e],[0 1e20],'-k')
% plot([-E_off_h-E_fermi_h, -E_off_h-E_fermi_h],[0 1e20],'-k')
% legend('Electrons (with broadening)','Holes (with broadening)','Electrons/Holes (no broadening)','Location','North')
% title('Convolution')
% xlim([-0.25, 2.5])
% ylim([0 1.6e19])
% hold off


