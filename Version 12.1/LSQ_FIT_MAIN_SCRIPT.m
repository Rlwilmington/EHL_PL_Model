%% MAIN SCRIPT - LSQ FIT OF EHL PL DATA

%be sure to have the following files also in the same folder:
%IBWread.m (necessary to import IGOR files)
%readIBWbinheader.m
%readIBWheaders.m
%D_e.m (weighted density of states convolved with Lorentzian degeneracy broadening)
%D_h.m
%PL_resid.m
%Ef_interpolation.m
%PLfitfun_int.m (electron-hole convolution to get PL)
%PL26Center_Corrected.ibw (data)
%EHL.mat (calibration)

%[f1, f2, f3, f4, f5] = StrainvTemp();

close all

[f1, f2, f3, f4, f5] = StrainvTemp_MAC();

%% OPTIONS

%filenum = 21; %pick which file to analyze (1-21)
cutoff_index = 0; %index left cutoff for data, zero for no cutoff
preview = 1; %if 1, input parameters will be plotted against raw data
useprev = 0; %if 1, previous fit parameters will be used as starting pars for this fit
maxevals = 2000; %max function evaluations allowed for LSQ fit
normalize = 0;

%Set Starting Parameters

UB = [6,  700,  0.8,  100, 100,   1, 0.05, 0.05]; %define fit bounds
LB = [0.3,  280,  0.4,  0,  0, 0, -0.01, -0.01];
deltaB = UB-LB;


par = LB + rand(1,length(UB)).*deltaB;

%par = [4.7687  603.8175    0.6588    8.9087    1.1942    0.5025   -0.0000    0.0022];

disp(par)

%% Import Raw Data

data = IBWread('PL26Center_Corrected_copy.ibw'); %raw intensity data vector
load('EHL.mat') %corresponding energies (calibration)

xdata = EHL(:,1);
ydata = data.y(:,filenum);


ydata = ydata*0.1; %just scale to 1, instead of 10


cutoff = find(xdata > cutoff_index); %set cutoff if necessary
xdata = xdata(cutoff:end);
ydata = ydata(cutoff:end); %scale data


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

if useprev == 1
    par = X;
    disp('Starting From Previous Endpoint')
end

const = [E_opt, E_bind, k_b, hbar, m0]; %constants necessary for fit model construction

%% Fit Preview

if preview == 1
    
    tic
    if normalize == 0
        [xfit_g, xfit_conv, yfit, u_vec, v_vec] = PLfitfun_int(par,const,f1,f2,f3,f4,f5);
    else
        [xfit_g, xfit_conv, yfit, u_vec, v_vec] = PLfitfun_int_norm(par,const,f1,f2,f3,f4,f5);
    end
    toc

    figure();
    hold on
    plot(xdata,ydata,'b*')
    plot(xfit_conv,yfit,'r-')
    legend('Raw Data', 'Trial Fit','Location','NorthWest')
    hold off
    
    figure();
    hold on
    plot(xfit_g,u_vec)
    plot(xfit_g,v_vec)
    plot(xfit_conv,yfit,'x')
    legend('D_e','D_h','Convolved Function')
    title('Convolution')
    hold off
    
    
    xfit_g_electrons = E_opt + xfit_g;
    xfit_g_holes = flipud(xfit_g);
    
    figure();
    hold on
    plot(xfit_g_electrons,u_vec)
    plot(xfit_g_holes,v_vec)
    legend('D_e','D_h','Convolved Function')
    title('Convolution')
    hold off
    
    drawnow
end

%% LSQ Nonlinear Fit
tic

minimizer = @(par) PL_resid(par,const,xdata,ydata,f1,f2,f3,f4,f5); %define single-var residual function to be minimized

fitops = optimoptions('lsqnonlin','maxfunevals',maxevals,'TolFun',1e-8,'TolX',1e-8,...
    'disp','iter-detailed','Algorithm','trust-region-reflective');

%options.Algorithm = 'sqp'; %'active-set'; %'interior-point'; %'levenberg-marquardt'; %select other algorithm

[X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(minimizer, par, LB, UB, fitops); %perform LSQ fit
toc
%% Results


res_y = PL_resid(X,const,xdata,ydata,f1,f2,f3,f4,f5);
[~, xfit_conv, yfit, ~, ~] = PLfitfun_int(X,const,f1,f2,f3,f4,f5);



% figure();
% hold on
% plot(xdata,ydata,'b.')
% plot(xfit_conv,yfit,'-r')
% legend('Raw Data', 'LSQ Fit')
% xlabel('Energy (eV)')
% ylabel('PL Intensity (Arb. Units)')
% hold off
% 
% figure();
% hold on
% plot(xdata,log(ydata),'b.')
% plot(xfit_conv,log(yfit),'-r')
% legend('Raw Data', 'LSQ Fit')
% xlabel('Energy (eV)')
% ylabel('Log PL Intensity (Arb. Units)')
% xlim([1.4 inf])
% ylim([-2 inf])
% hold off
% 
% figure();
% hold on
% title('LSQ Fit Residuals')
% plot(xdata,res_y,'*');
% plot([xdata(1), xdata(end)],[0 0],'k-')
% hold off

J = JACOBIAN;
J = full(J);
dof = length(xdata) - length(X);
mse = (1/dof)*sum(res_y.^2); %variance of residuals (assuming ~uniform)
cov = (J.'*J)*mse;
SE = sqrt(diag(cov));
CI = nlparci(X,res_y,'covar',cov);

disp('10^13 Charge Carrier Density per Volume:')
disp(X(1))

disp('Temperature (K):')
disp(X(2))

disp('Electrons: \Gamma = alpha(E-Ef)^2 + C')
disp('alpha')
disp(X(4))
disp('C')
disp(X(7))

disp('Holes: \Gamma = beta(E-Ef)^2 + D')
disp('beta')
disp(X(5))
disp('D')
disp(X(8))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if X(2) > 294
    strain = polyval(f1,X(2));
else
    strain = polyval(f1,294);
end


disp('Sample Strain (Lattice Expansion)')
disp(strain)
E_off_h = polyval(f3,strain); %eV
m_e_K = polyval(f5,strain)*m0; %effective masses (highest fluence)
m_h_K = polyval(f4,strain)*m0;
m_h_G = polyval(f2,strain)*m0;
g_e_K = (degen_e_K*m_e_K)/(pi*hbar^2); %density of states
g_h_K = (degen_h_K*m_h_K)/(pi*hbar^2);
g_h_G = (degen_h_G*m_h_G)/(pi*hbar^2);

const = [E_opt, E_bind, k_b, hbar, m0]; %constants necessary for fit model construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Electron and Hole Fermi Level (eV)')
E_fermi_e = F1(X(2),X(1).*1e13,const,g_e_K);
E_fermi_h = F2(X(2),X(1).*1e13,const,g_h_K,g_h_G,E_off_h);
display(E_fermi_e)
display(E_fermi_h)

a = 1e-4; %units
T = X(2);

n_h_K = a*g_h_K*k_b*T*log(exp(E_fermi_h/(k_b*T)) + 1);

n_h_G = a*g_h_G*k_b*T*log(exp((E_fermi_h + E_off_h)/(k_b*T)) + 1);

disp('Total Hole Density')
disp(X(1)*1e13)
disp(n_h_K + n_h_G)

disp('K-Valley Hole Density')
disp(n_h_K)

disp('G-Valley Hole Density')
disp(n_h_G)

disp('Bandgap Renormalization Energy (eV):')
disp(X(3))

disp('Predicted BGR')
a_bohr = 6.02e-8;
bgr = @(n) 3.1*(((n*1e13).*a_bohr^2).^(1/3)).*(E_bind);
disp(bgr(X(1)))

KKgap = @(T) (353e-6)*(T-294); %from Avinash calculations (E_thermal)


if T > 294
    E_thermal = KKgap(T);
else
    E_thermal = 0;
end

disp('Thermal Bandgap Reduction')
disp(E_thermal)


