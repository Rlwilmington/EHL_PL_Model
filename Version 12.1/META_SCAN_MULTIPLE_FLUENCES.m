%% MULTIPLE FLUENCE SCANNING DATA

%% OPTIONS

maxevals = 1500; %max function evaluations allowed for LSQ fit
files = 21;
scans = 5;
flux = 0.05; %how much to vary input par per scan
skip = 1;

%% Import Raw Data

data = IBWread('PL26Center_Corrected_copy.ibw'); %raw intensity data vector
load('EHL.mat') %corresponding energies (calibration)

xdata = EHL(:,1);

xdata = xdata(1:skip:end);

[f1, f2, f3, f4, f5] = StrainvTemp();


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
E_bind = 0.42; %eV, note that exciton binding energy estimates vary 0.4-0.8eV (0.53 suggested by Robert Younts, 0.42 by Bataller et al)

const = [E_opt, E_bind, k_b, hbar, m0]; %constants necessary for fit model construction



%% Multiple Fluence Scans

allscans = zeros(files,8,scans);
k = 1;
UB = [6,  700,  0.9,  30, 0.2,   10]; %define fit bounds
LB = [0.1,  200,  0.4,  0,  0.001, 0];
deltaB = (UB-LB);

for i = 1:files
    
    filenum = i;
    j = 1;
    ydata = data.y(:,filenum);
    ydata = 0.1*ydata(1:skip:end);
    
    filedatamatrix = zeros(scans,8);
    
    while j <= scans
        
        disp([i, j])
        %par = normrnd(par0,range); %randomize par for each scan
        par = LB + rand(1,6).*deltaB;
        
        try
            [X,RESNORM] = PL_fit_generator(par,const,xdata,ydata,f1,f2,f3,f4,f5,maxevals);
            newscan = horzcat(filenum,X,RESNORM);
            filedatamatrix(j,:) = newscan;
            j = j + 1;
        catch
            disp('Scan Failed, Trying Again...')
        end
        
    end
    
    allscans(i,:,:) = filedatamatrix.';
    
end


%%

allscans2D = zeros(files*scans,8);
ind = 1:files:scans*files;

for i = 1:scans
    k = ind(i);
    allscans2D(k:k+files-1,:) = allscans(:,:,i);
end

allscans2D = sortrows(allscans2D,1);

%%

Y = allscans2D;

% allscans2D = allscans2D(allscans2D(:,3) > 301,:);

Y = zeros(files,8);
for i = 1:files
    scan_list = allscans2D(allscans2D(:,1) == i,:);
    scan_list = sortrows(scan_list,8);
    try
        Y(i,:) = scan_list(1,:);
    catch
    end
end


%% BGR check

x = Y(:,2);
z = Y(:,4);

a_bohr = 6.5e-8;
fix = 0;
bgr = @(n) 2.7*(((n*1e13).*a_bohr^2).^(1/3)).*(E_bind+fix) - fix;
nspace = linspace(0.1,5,1000);
bgrvals = bgr(nspace);

% min_g = @(x,par) par(1)*(((x*1e13).*par(2)).^(1/3)).*(0.42+par(3)) - par(3);
% 
% par = [3.1,6.5e-8,0.01];
% UB = [5,10e-8,0.3]; %define fit bounds
% LB = [1,3e-8,-0.3];
% 
% fitops = optimoptions('lsqnonlin','maxfunevals',maxevals,'TolFun',1e-12,'TolX',1e-12); %,'Display','iter-detailed','Algorithm','trust-region-reflective');
% 
% [X, RESNORM] = lsqnonlin(min_g, par, LB, UB, fitops);

figure();
hold on
plot(nspace,bgrvals,'r')
plot(x,z,'*')
hold off





%% META OUTPUT

E_off = [
    564.453125, -0.05789473684210525
    818.359375, -0.050375939849624074
    1056.6406249999998, -0.054385964912280704
    1306.6406249999998, -0.051629072681704274
    1541.0156249999998, -0.04160401002506267
    1791.0156249999998, -0.013533834586466176
    2037.1093749999993, -0.010526315789473675
    2275.390624999999, -0.0027568922305764437
    2517.578125, 0.01177944862155389
    2763.671874999999, 0.026315789473684216
    3009.765624999999, 0.03834586466165414
    3251.953124999999, 0.03007518796992481
    3501.953124999999, 0.047368421052631574
    3736.328124999999, 0.047368421052631574
    3986.328124999999, 0.057393483709273184];

list = 1:15;

fluence_vec = interp1(list,1e-3*E_off(:,1),Y(:,1),'linear','extrap');

%%

figure();
semilogy(fluence_vec,Y(:,8),'*')
title('Fit Quality')
ylabel('MSE')
xlabel('Power Density (kW/cm^2)')

figure();
plot(fluence_vec,Y(:,2),'*')
title('Number Density vs Fluence')
xlabel('Power Density (kW/cm^2)')
ylabel('Number Density (1e13/cm^2)')

figure();
plot(fluence_vec,Y(:,3),'*')
title('Temperature vs Fluence')
xlabel('Power Density (kW/cm^2)')
ylabel('Temperature (K)')

figure();
plot(fluence_vec,Y(:,4),'*')
title('BGR vs Fluence')
xlabel('Power Density (kW/cm^2)')
ylabel('Energy (eV)')

figure();
plot(fluence_vec,Y(:,5),'*')
title('Alpha vs Fluence')
xlabel('Power Density (kW/cm^2)')
ylabel('Alpha (arb.)')

figure();
plot(fluence_vec,Y(:,6),'*')
title('Beta vs Fluence')
xlabel('Power Density (kW/cm^2)')
ylabel('Beta (arb.)')

figure();
plot(fluence_vec,Y(:,7),'*')
title('Quasi-Quantum Efficiency vs Fluence')
xlabel('Power Density (kW/cm^2)')
ylabel('Quasi-Quantum Efficiency (arb.)')

figure();
plot(fluence_vec,Y(:,8),'*')
title('C vs Fluence')
xlabel('Power Density (kW/cm^2)')
ylabel('C (arb.)')

figure();
plot(fluence_vec,Y(:,9),'*')
title('D vs Fluence')
xlabel('Power Density (kW/cm^2)')
ylabel('D (arb.)')


% figure();
% plot(Y(:,3),Y(:,8),'*')
% title('Temperature vs Fit Quality')
% 
% figure();
% plot(Y(:,2),Y(:,8),'*')
% title('Number Density vs Fit Quality')

% for i = 1:files
%     figure();
%     plotdata = Y(Y(:,1) == i,:);
%     plot(plotdata(:,2),plotdata(:,8),'*')
%     ttl = strcat('Number Density vs Fit Quality: ', num2str(i));
%     title(ttl)
% end

%% Fermi Levels, E_Thermal


Z = zeros(length(Y),5);

const = [E_opt, E_bind, k_b, hbar, m0]; %constants necessary for fit model construction

for i = 1:length(Y)
    
    n = Y(i,2);
    T = Y(i,3);
    
    strain = f1(T);
    E_off_h = f3(strain); %eV
    m_e_K = f5(strain)*m0; %effective masses (highest fluence)
    m_h_K = f4(strain)*m0;
    m_h_G = f2(strain)*m0; %1.9*m0;
    g_e_K = (degen_e_K*m_e_K)/(pi*hbar^2); %density of states
    g_h_K = (degen_h_K*m_h_K)/(pi*hbar^2);
    g_h_G = (degen_h_G*m_h_G)/(pi*hbar^2);

    E_fermi_e = F1(T,n.*1e13,const,g_e_K);
    E_fermi_h = F2(T,n.*1e13,const,g_h_K,g_h_G,E_off_h);

    a = 1e-4; %units

    n_h_K = a*g_h_K*k_b*T*log(exp(E_fermi_h/(k_b*T)) + 1);

    n_h_G = a*g_h_G*k_b*T*log(exp((E_fermi_h + E_off_h)/(k_b*T)) + 1);

    KKgap = @(T) (-353e-6)*T + 2.0118; %from Avinash calculations (E_thermal)
    E_thermal = E_opt - KKgap(T);

    Z(i,:) = [E_fermi_e, E_fermi_h, n_h_K*1e-13, n_h_G*1e-13, E_thermal];
end

%% Plots

figure();
hold on
plot(fluence_vec,Z(:,1),'r*')
plot(fluence_vec,Z(:,2),'b*')
title('Fermi Level vs Fluence')
legend('Electron Ef','Hole Ef')
xlabel('Power Density (kW/cm^2)')
ylabel('Energy (eV)')
hold off

figure();
hold on
plot(fluence_vec,Z(:,3),'r*')
plot(fluence_vec,Z(:,4),'b*')
plot(fluence_vec,Y(:,2),'k*')
xlabel('Power Density (kW/cm^2)')
ylabel('Number Density (1e13/cm^2)')
title('Number Density per Valley vs Fluence')
legend('K-Valley','\Gamma Valley', 'Total')
hold off


%% Histograms

for i = 1:files
    
    histdata = Y(Y(:,1) == i,:);
    hist(histdata(:,2))
    title(num2str(i))
    drawnow
    pause
    
end
