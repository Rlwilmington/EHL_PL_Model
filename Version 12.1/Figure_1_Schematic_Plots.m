

clear fit M1 M2

fontsize = 16;

[f1, f2, f3, f4, f5] = StrainvTemp_MAC();

m0_oldunits = 1e6*0.510998950; %eV/c^2
c = 299792458; %m/s
m0 = (1/c^2)*m0_oldunits;
hbar = 6.582119569e-16; %eV*s
k_b = 8.617333262145e-5; %eV/K

degen_e_K = 2; %band degeneracy
degen_h_K = 2;
degen_h_G = 1;

E_opt = 1.89; %eV
E_bind = 0.421; %eV, note that exciton binding energy estimates vary 0.4-0.8eV (0.53 suggested by Robert Younts, 0.42 by Bataller et al)

startfile = 1;

% E_off = [
%     564.453125, -0.05789473684210525
%     818.359375, -0.050375939849624074
%     1056.6406249999998, -0.054385964912280704
%     1306.6406249999998, -0.051629072681704274
%     1541.0156249999998, -0.04160401002506267
%     1791.0156249999998, -0.013533834586466176
%     2037.1093749999993, -0.010526315789473675
%     2275.390624999999, -0.0027568922305764437
%     2517.578125, 0.01177944862155389
%     2763.671874999999, 0.026315789473684216
%     3009.765624999999, 0.03834586466165414
%     3251.953124999999, 0.03007518796992481
%     3501.953124999999, 0.047368421052631574
%     3736.328124999999, 0.047368421052631574
%     3986.328124999999, 0.057393483709273184];
% 
% file_list = 1:15;

%power_vec = (1e-3*interp1(file_list,E_off(:,1),1:21,'linear','extrap').');
power_vec = 1e-3*(570.302+(0:20).*243.621).';

%%
allentries = zeros(21,18);
Mxdata = zeros(21,1109);
Mydata = zeros(21,1109);


for i = startfile:21
    filenum = i;
    filestring = num2str(filenum);
    filename = strcat('montecarloerrorsaveVersion12.11_April1720_file',filestring,'.mat');

    S = load(filename);
    
    entry = [filenum,S.X,S.RESNORM,S.s_mat.'];
    allentries(i,:) = entry;
    Mxdata(i,:) = S.xdata;
    Mydata(i,:) = S.ydata;

end

% for i = startfile:21 %makes waterfall plot
%     Mydata(i,:) = Mydata(i,:) + 0.25*(i-startfile);
% end

%%


% Mxfit = Mxfit(:,fitplotbounds(1):fitplotbounds(2));
% Myfit = Myfit(:,fitplotbounds(1):fitplotbounds(2));


start_vec = zeros(21,1);

%%




%% Load Data

load LatticeExpVSPower.mat LatticeExpVSPower 
load TempVsPower.mat TempVsPower

M1 = sortrows(LatticeExpVSPower,1);
M2 = sortrows(TempVsPower,1);

var1 = M1(:,1);
var2 = M2(:,1);

fit = fit(var1,var2,'poly1');
f = [fit.p1, fit.p2];

%% Calibrate Sample Temperature to Percent Strain Using Raman Data

LatticeExpVSPower(:,1) = flipud(LatticeExpVSPower(:,1));
tempStrainCoeff = polyfit(TempVsPower(1:5,1),LatticeExpVSPower(1:5,1),1);

%%Discretize Space
numPts = 10;
T = linspace(295,525,numPts);
latticeExp = polyval(tempStrainCoeff,T);

%% Determine Energy Shifts of K-CB, Q-CB, G-VB relative to the K-VB for the discretized T/Lattice Expansion
% Code Written by Dr. Robert Younts
load EnergyOffsetStrain.mat energyCoeff
%energyOffset = table;
energyOffset = table('Size',[numPts 4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'Strain','K_CB','Q_CB','G_VB'});
energyOffset.Strain = latticeExp';
energyOffset.K_CB = polyval(energyCoeff(:,1),energyOffset.Strain);
energyOffset.Q_CB = polyval(energyCoeff(:,2),energyOffset.Strain);
energyOffset.G_VB = polyval(energyCoeff(:,3),energyOffset.Strain);

M = table2array(energyOffset);

%%

temps = polyval(f,M(:,1));

KKgap = @(T) (-353e-6)*(T-293); % + 2.01; %1.9934;

KKenergyOffset = E_opt - KKgap(temps);



%% PLOTS

low = M(1,1);
high = M(end,1);
pad = 0.06;




figure('Renderer', 'painters', 'Position', [10 10 1400 450])
%figure('units','inch','position',[0,0,3.3,2*3.3/3]);

axes('Position',[0.001 0.001 0.333 0.998],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0.333 0.001 0.333 0.998],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0.666 0.001 0.333 0.998],'xtick',[],'ytick',[],'box','on','handlevisibility','off')

left_color = [0 0 0];
right_color = [0 0 0];
ax_height_start = 0.18;
ax_height = 0.77;

h1 = axes('Position',[0.05 ax_height_start 0.255 ax_height]);
set(gca,'FontSize',fontsize,'CLim',[power_vec(startfile) power_vec(21)])

m_color = (256-1)/(power_vec(21)-power_vec(1));
b_color = 1 - m_color*power_vec(1);

color_list = round(m_color.*power_vec + b_color);


colormap(jet);
A = colormap;
custcolormap = A(color_list(startfile:21),:);

ColorOdrDef = get(gca,'ColorOrder'); %7x3 RGB array
ColorOdrCustom = flipud(custcolormap);
%ColorOdrCustom = custcolormap;
set(gca,'ColorOrder',ColorOdrCustom);

hold on
for i = fliplr(startfile:21)
%for i = startfile:2:17
%for i = startfile:21
    filenum = i;
    p1(i) = plot(h1,Mxdata(filenum,:),Mydata(filenum,:),'o','MarkerSize',2);
end

xlim([1.45 2.05])
ylim([-0.05 1.05])
ylabel('PL Intensity (arb.)')
xlabel('Energy (eV)')
c = colorbar(gca);
c.Label.String = 'Power Density (kW/cm^2)';
hold off


h3 = axes('Position',[0.72 ax_height_start 0.26 ax_height]);
set(gca,'FontSize',fontsize)
hold on
plot(h3,M(:,1),M(:,2),'-o','LineWidth',2);
plot(h3,M(:,1),M(:,3),'-^','LineWidth',2);
plot(h3,M(:,1),M(:,2) - M(:,4),'-s','LineWidth',2);
plot(h3,M(:,1),M(:,3) - M(:,4),'-d','LineWidth',2);
ylabel('Energy Gap (eV)')
xlabel('Sample Strain (%)')
legend('K(VB)-K(CB)','K(VB)-Q(CB)','\Gamma(VB)-K(CB)','\Gamma(VB)-Q(CB)','Location','SouthWest')
xlim([low-pad high+pad])
ylim([1.35, 2.0])
dim = [0.007 0.88 0.1 0.1];
str = {'(a)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.340 0.88 0.1 0.1];
str = {'(b)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.673 0.88 0.1 0.1];
str = {'(c)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
hold off


set(gcf,'defaultAxesColorOrder',[left_color; right_color]);

h2 = axes('Position',[0.40 ax_height_start+0.01 0.20 ax_height]);
set(gca,'FontSize',fontsize)
range = 9:15;
ylim_temp = [250, 650];
ylim_strain = polyval(f1,ylim_temp);
strain_list = polyval(f1,allentries(range,3));
yyaxis left
plot(power_vec(range),allentries(range,3),'-bo','LineWidth',2);
ylim(ylim_temp)
xlim([2.3, 4.2])
ylabel('Temperature (K)')

yyaxis right
plot(power_vec(range),strain_list,'-bo','LineWidth',1);
ylim(ylim_strain)
ylabel('Strain (%)')
xlabel('Power Density (kW/cm^2)')


saveas(gcf,'fig_fig1.png')

%%



