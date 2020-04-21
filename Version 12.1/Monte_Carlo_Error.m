%% MONTE CARLO ERROR ESTIMATION FOR LSQ FIT RESULT

%% OPTIONS

%filenum = 19; %pick which file to analyze (1-21)
cutoff_index = 0; %index left cutoff for data, zero for no cutoff
find_fit = 0;
skip = 1; %only use every __ data points to speed up fit
maxevals = 700; %max function evaluations allowed for LSQ fit
trials = 250; %false datasets created for error estimation
bins = 60; %histogram bins for parameter error estimates
start = 1;

%Set Starting Parameters
UB = [6,  700,  0.8,  100, 100,   1, 0.05, 0.05]; %define fit bounds
LB = [0.3,  280,  0.4,  0,  0, 0, -0.01, -0.01];
deltaB = UB-LB;

[f1, f2, f3, f4, f5] = StrainvTemp_MAC();
%% Import Raw Data

data = IBWread('PL26Center_Corrected_copy.ibw'); %raw intensity data vector
load('EHL.mat') %corresponding energies (calibration)

xdata = EHL(:,1);
ydata = data.y(:,filenum);
ydata = 0.1*ydata; %normalize


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

%% Find Fit

if find_fit == 1
    minimizer = @(par) PL_resid(par,const,xdata,ydata,F1,F2);
    %par0 = [A,n,T,E_BGR,alpha,beta]; %init

    fitops = optimoptions('lsqnonlin','maxfunevals',maxevals,'TolFun',1e-8,'TolX',1e-8,...
        'Display','iter-detailed','Algorithm','trust-region-reflective');
    %options.Algorithm = 'levenberg-marquardt'; %select other algorithm
    X = lsqnonlin(minimizer, par, LB, UB, fitops);
end

[~, xfit_conv, yfit, ~, ~] = PLfitfun_int(X,const,f1,f2,f3,f4,f5);
res_y = PL_resid(X,const,xdata,ydata,f1,f2,f3,f4,f5);
yfit_interp = interp1(xfit_conv,yfit,xdata);

% figure();
% hold on
% plot(xdata,ydata)
% plot(xfit_conv,yfit)
% hold off

dof = length(xdata) - length(X);
mse = (1/dof)*sum(res_y.^2); %variance of residuals (assuming ~uniform)
rmse = sqrt(mse);

%% False Data Creation

X_err = zeros(trials,length(X)+1);
par = X;
sigma = abs(par)*0.001;

tic
for i = start:trials    
    
    par_trial = LB + rand(1,length(UB)).*deltaB;
    err_data = normrnd(yfit_interp,rmse); %generate artificial data
    try
        [X_trial,RESNORM] = PL_fit_generator(par_trial,const,xdata,err_data,f1,f2,f3,f4,f5,maxevals);
        X_err(i,:) = [X_trial,RESNORM];
    catch
        X_err(i,:) = ones(1,9);
    end
    disp(i)
end
toc

%% Fit Parameter Histograms to Gaussians

X_err_clean = X_err(X_err(:,9) < 1,:);
s_mat = zeros(length(X),1);

for i = 1:length(X)
    
    fitbins = bins; %length(counts) - length(find(counts == 0));

    fitdata = X_err_clean(:,i);
    plotmin = min(fitdata);
    plotmax = max(fitdata);
    counts = histcounts(fitdata,fitbins,'BinLimits',[plotmin,plotmax]);
    norm = (1/max(counts))*counts;

    fd = fitdist(fitdata,'Normal');

%     figure()
%     fit = histfit(fitdata,fitbins,'Normal');
%     legend('Data','Fit')
%     %xlim([plotmin plotmax])
%     xlabel('x')
%     ylabel('Counts')
    
    s_mat(i) = fd.sigma;
end

%% Find Confidence Intervals

CI = zeros(length(X),2);

for i = 1:length(X)
    
    CI(i,:) = [X(i) - s_mat(i), X(i) + s_mat(i)];
    
end
