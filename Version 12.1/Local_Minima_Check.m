%% Evaluate Quality of Local Minima


%% Options

filenum = 21; %pick which file to analyze (1-21)
cutoff_index = 0; %index left cutoff for data, zero for no cutoff
preview = 0; %if 1, input parameters will be plotted against raw data
useprev = 0; %if 1, previous fit parameters will be used as starting pars for this fit
skip = 1; %only use every __ data points to speed up fit
maxevals = 300; %max function evaluations allowed for LSQ fit

%Set Starting Parameters
par0 = [4.0,  440,    0.69,    2.5,    0.05];
range = par0*0.1;


% par0 = [0.0375,   2.3,  405,    0.70,    1.7,    0.051]; %file 21 good fit, %par0 = [A,n,T,E_BGR,alpha,beta];

%% Import Raw Data

data = IBWread('PL26Center_Corrected_copy.ibw'); %raw intensity data vector
load('EHL.mat') %corresponding energies (calibration)

xdata = EHL(:,1);
ydata = data.y(:,filenum);
[maxval,maxind] = max(ydata);
scale = mean(ydata(maxind-15:maxind+15));
ydata = (1/scale)*ydata; %normalize

cutoff = find(xdata > cutoff_index); %set cutoff if necessary
xdata = xdata(cutoff:end);
ydata = ydata(cutoff:end);

%Edit Data (?)

xdata = vertcat(xdata(1:865),xdata(877:end));
ydata = vertcat(ydata(1:865),ydata(877:end));

%% Set Up Model Fit

degen_e_K = 2; %band degeneracy
degen_h_K = 2;
degen_h_G = 1;
m0_oldunits = 1e6*0.510998950; %eV/c^2
c = 299792458; %m/s
m0 = (1/c^2)*m0_oldunits;
m_e_K = 0.44*m0; %effective masses (highest fluence)
m_h_K = 0.525*m0;
m_h_G = 1.9*m0;
hbar = 6.582119569e-16; %eV*s
k_b = 8.617333262145e-5; %eV/K
E_off_h = 0.055; %eV

g_e_K = (degen_e_K*m_e_K)/(pi*hbar^2); %density of states
g_h_K = (degen_h_K*m_h_K)/(pi*hbar^2);
g_h_G = (degen_h_G*m_h_G)/(pi*hbar^2);

N = length(xdata); %should be 1109 points

x0 = xdata(1);
xf = xdata(N);
xrange = xf-x0;
E_opt = 1.89; %eV
E_bind = 0.42; %eV, note that exciton binding energy estimates vary 0.4-0.8eV (0.53 suggested by Robert Younts, 0.42 by Bataller et al)

%Set Starting Parameters

%par0 = [1.6,   1.34,  390,    0.28,    0.007,    0.05];
%par0 = [1.5197,1.3425,390.6286,0.2707,0.0202,0.0590]; %file 20
%par0 = [2,   1.15,  450,    0.28,    0.08,    0.02];
%par0 = [1.8926,    1.1306,  445.0196,    0.2752,    0.0984,    0.0543]; %interesting fit file 20
%par0 = [0.4,  0.95,  350,    0.18,    0.4,    0.055];

par = par0;

if useprev == 1
    par = X;
    disp('Starting From Previous Endpoint')
end

const = [E_opt, E_bind, k_b, g_e_K, g_h_K, g_h_G, hbar, E_off_h]; %constants necessary for fit model construction

xdata_reduc = xdata(1:skip:N);
ydata_reduc = ydata(1:skip:N);

num = length(xdata_reduc);

[F1, F2] = Ef_interpolation(const); %build fermi level interpolator functions, ex: electron fermi level = F1(temp, number density)




%% Multistart

nsteps = 100; %try 100 different starting number densities
trials = 10; %try 10 trials at each n
fitdata = zeros(nsteps*trials, length(par0) + 1); %record all data in one giant matrix
k = 1;

for j = 1:nsteps
    
    par0(1) = 1.0 + j*(1/nsteps); %n will range from 1 to 2

    for i = 1:trials

        par = [par0(1), normrnd(par0(2:end),range(2:end))]; %create randomized start point

        try
            [X, RESNORM] = PL_fit_generator(par,const,xdata_reduc,ydata_reduc,F1,F2,maxevals);
        catch
            X = zeros(1,length(par0));
            RESNORM = 0;
        end

        fitdata(k,:) = horzcat(X, RESNORM);
        disp(fitdata(k,:))
        k = k+1;

    end

end

%% Remove blank rows


fitdata_sort = sortrows(fitdata);
fitdata_sort( all(~fitdata_sort,2), : ) = []; %remove empty rows

bad_temp = find(fitdata_sort(:,2) < 400);
fitdata_sort(bad_temp,:) = [];

figure();
semilogy(fitdata_sort(:,1),fitdata_sort(:,6),'*')
xlabel('n')
ylabel('MSE')

figure();
plot(fitdata_sort(:,1),fitdata_sort(:,6),'*')
xlabel('n')
ylabel('MSE')

%% Bin Data

bins = 100;
bin_min = 0.7;
limit = 2.3;
step = (limit-bin_min)/bins;
bin_max = bin_min + step;
x = zeros(bins,1);
y = zeros(bins,1);
z = zeros(bins,1);

for i = 1:bins
    
    bin_ind = find(fitdata_sort(:,1) >= bin_min & fitdata_sort(:,1) < bin_max);
    x(i) = (bin_min+bin_max)/2;
    val = min(fitdata_sort(bin_ind,6));
    if isempty(val)
        y(i) = -1;
        z(i) = -1;
    else
        y(i) = val;
        bgr_ind = find(fitdata(:,6) == val);
        z(i) = fitdata(bgr_ind,3);
    end
    
    bin_min = bin_min + step;
    bin_max = bin_max + step;
end

nodata = find(y == -1);

x(nodata) = [];
y(nodata) = [];
z(nodata) = [];

plot(x,y,'x')

%% BGR check

a_bohr = 6.5e-8;

bgr = @(n) 3.1*(((n*1e13).*a_bohr^2).^(1/3)).*(E_bind);

nspace = linspace(0.7,3,1000);
bgrvals = bgr(nspace);


hold on
plot(nspace,bgrvals,'r')
plot(x,z - E_bind,'*')
hold off


