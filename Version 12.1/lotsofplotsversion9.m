%%

%MOST RECENT PLOTTER, INCLUDES ERRORBARS AND ROBERT CHANGES
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

startfile = 8;

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

file_list = 1:15;


%%
allentries = zeros(21,18);
Mxdata = zeros(21,1109);
Mydata = zeros(21,1109);
Mxfit = zeros(21,401);
Myfit = zeros(21,401);



for i = startfile:21
    filenum = i;
    filestring = num2str(filenum);
    filename = strcat('montecarloerrorsaveVersion12.11_April1720_file',filestring,'.mat');

    S = load(filename);
    
    entry = [filenum,S.X,S.RESNORM,S.s_mat.'];
    allentries(i,:) = entry;
    Mxdata(i,:) = S.xdata;
    Mydata(i,:) = S.ydata;
    Mxfit(i,:) = S.xfit_conv;
    Myfit(i,:) = S.yfit;

end

peakheight = zeros(21,1);

for i = startfile:21
    
    peakheight(i) = max(Mydata(i,:));
    
end

peakheight = (1/max(peakheight))*peakheight;

for i = startfile:21
    Myfit(i,:) = Myfit(i,:) + 0.25*(i-startfile);
    Mydata(i,:) = Mydata(i,:) + 0.25*(i-startfile);
end

%power_vec = (1e-3*interp1(file_list,E_off(:,1),1:21,'linear','extrap').');
power_vec = 1e-3*(570.302+(0:20).*243.621).';

%%

N = 1109;
n = 500;

T = readtable('alldataVersion10Feb1820.txt');
alldata = table2array(T);
Mxfit10 = alldata(2*N+1:2*N+n,:);
Myfit10 = alldata(2*N+n+1:2*(N+n),:);
fitdata = alldata(end-5:end,:);

start_vec = zeros(21,1);

for i = startfile:21
    ind_list = find(Myfit10(:,1) == 0);
    start_vec(i) = ind_list(end);
    Myfit10(:,i) = Myfit10(:,i) + 0.25*(i-startfile);
end

%% OLD DUAL WATERFALL PLOT

% line_val = 1.3;
% 
% 
% figure('Renderer', 'painters', 'Position', [10 10 1200 600])
% axes('Position',[0.001 0.001 0.498 0.998],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
% axes('Position',[0.5 0.001 0.498 0.998],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
% 
% h1 = axes('Position',[0.07 0.13 0.38 0.8]);
% set(gca,'FontSize',fontsize,'CLim',[power_vec(startfile) power_vec(21)])
% 
% m_color = (256-1)/(power_vec(21)-power_vec(1));
% b_color = 1 - m_color*power_vec(1);
% 
% color_list = round(m_color.*power_vec + b_color);
% 
% 
% colormap(jet);
% A = colormap;
% custcolormap = A(color_list(startfile:21),:);
% 
% ColorOdrDef = get(gca,'ColorOrder'); %7x3 RGB array
% ColorOdrCustom = custcolormap;
% set(gca,'ColorOrder',ColorOdrCustom);
% 
% 
% hold on
% for i = startfile:21
%     filenum = i;
%     p1(i) = plot(Mxdata(filenum,:),Mydata(filenum,:),'o','MarkerSize',2);
% end
%     
%     
% for i = startfile:21
%     filenum = i;
%     p2(i) = plot(Mxfit(filenum,:),Myfit(filenum,:),'k:','LineWidth',line_val);
% end
% 
% 
% %legend([p2(21), p2(20), p2(19), p2(18), p2(17), p2(16), p2(15), p2(14), p2(13), p2(12), p2(11), p2(10), p2(9), p2(8), p2(7), p2(6), p2(5), p2(4)],leg)
% xlim([1.4 2.1])
% ylim([-0.25 4.4])
% ylabel('PL Intensity (arb.)')
% xlabel('Energy (eV)')
% c = colorbar(gca);
% c.Label.String = 'Power Density (kW/cm^2)';
% legend([p2(startfile)],'Model Fit')
% hold off
% 
% h2 = axes('Position',[0.57 0.13 0.38 0.8]);
% set(gca,'FontSize',fontsize,'CLim',[power_vec(startfile) power_vec(21)])
% set(gca,'ColorOrder',ColorOdrCustom);
% 
% hold on
% for i = startfile:21
%     filenum = i;
%     p1(i) = plot(Mxdata(filenum,:),Mydata(filenum,:),'o','MarkerSize',2);
% end
% 
% for i = startfile:21
%     filenum = i;
%     p2(i) = plot(Mxfit10(start_vec(i):end,filenum),Myfit10(start_vec(i):end,filenum),':k','LineWidth',line_val);
% end
% 
% xlim([1.4 2.1])
% ylim([-0.25 4.4])
% ylabel('PL Intensity (arb.)')
% xlabel('Energy (eV)')
% legend([p2(startfile)],'Model Fit')
% c = colorbar(gca);
% c.Label.String = 'Power Density (kW/cm^2)';
% dim = [0.007 0.88 0.1 0.1];
% str = {'(a)'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
% dim = [0.507 0.88 0.1 0.1];
% str = {'(b)'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
% hold off
% 
% 
% saveas(gcf,'fig_fitwith&withoutL.png')

%%

line_val = 1.3;


figure('Renderer', 'painters', 'Position', [10 10 900 800])
axes('Position',[0.001 0.001 0.667 0.999],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0.668 0.666 0.333 0.334],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0.668 0.333 0.333 0.333],'xtick',[],'ytick',[],'box','on','handlevisibility','off')

h1 = axes('Position',[0.09 0.1 0.51 0.82]);
%set(gca,'FontSize',fontsize,'CLim',[power_vec(startfile) power_vec(21)])
color_power_start = 1.5;
color_power_end = 5.75;
set(gca,'FontSize',fontsize,'CLim',[color_power_start color_power_end])

m_color = (256-1)/(color_power_end-color_power_start);
b_color = 1 - m_color*color_power_start;

color_list = round(m_color.*power_vec + b_color);


colormap(jet);
A = colormap;
custcolormap = A(color_list(startfile:21),:);

ColorOdrDef = get(gca,'ColorOrder'); %7x3 RGB array
ColorOdrCustom = custcolormap;
set(gca,'ColorOrder',ColorOdrCustom);


hold on
for i = startfile:21
    filenum = i;
    p1(i) = plot(Mxdata(filenum,:),Mydata(filenum,:),'o','MarkerSize',3);
end
    
    
for i = startfile:21
    filenum = i;
    p2(i) = plot(Mxfit(filenum,:),Myfit(filenum,:),'k:','LineWidth',line_val);
end


%legend([p2(21), p2(20), p2(19), p2(18), p2(17), p2(16), p2(15), p2(14), p2(13), p2(12), p2(11), p2(10), p2(9), p2(8), p2(7), p2(6), p2(5), p2(4)],leg)
xlim([1.4 2.1])
ylim([-0.25 4.4])
ylabel('PL Intensity (arb.)')
xlabel('Energy (eV)')
c = colorbar(gca);
c.Label.String = 'Power Density (kW/cm^2)';
legend([p2(startfile)],'Model Fit','FontSize',12)
hold off

h2 = axes('Position',[0.73 0.727 0.26 0.26]);
set(gca,'FontSize',10,'CLim',[power_vec(startfile) power_vec(21)])

hold on
filenum = 21;
p1 = plot(Mxdata(filenum,:),Mydata(filenum,:)-3.25,'.','MarkerSize',10);
p2 = plot(Mxfit10(start_vec(i):end,filenum),Myfit10(start_vec(i):end,filenum)-3.25,'r--','LineWidth',1.5); 
%3.25 to bring scale back to zero, only showing 1 plot instead of waterfall

xlim([1.44 2.07])
ylim([0 0.95])
ylabel('PL Intensity (arb.)')
xlabel('Energy (eV)')
legend('Raw Data','Simplified Fit','Location','NorthEast','FontSize',8)
dim = [0.002 0.905 0.1 0.1];
str = {'(a)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.666 0.905 0.1 0.1];
str = {'(b)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.666 0.571 0.1 0.1];
str = {'(c)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
hold off

h3 = axes('Position',[0.73 0.395 0.26 0.26]);
set(gca,'FontSize',10,'CLim',[power_vec(startfile) power_vec(21)])

hold on
filenum = 21;
p1 = plot(Mxdata(filenum,:),Mydata(filenum,:)-3.25,'.','MarkerSize',10);
p3 = plot(Mxfit(filenum,:),Myfit(filenum,:)-3.25,'r--','LineWidth',1.5);

xlim([1.44 2.07])
ylim([0 0.95])
ylabel('PL Intensity (arb.)')
xlabel('Energy (eV)')
legend('Raw Data','Complete Fit','Location','NorthEast','FontSize',8)
hold off

saveas(gcf,'fig_fig3','svg')
saveas(gcf,'fig_fig3','png')



%%

fontsize = 16;

figure('Renderer', 'painters', 'Position', [10 10 650 800])
axes('Position',[0.001 0.667 0.998 0.333],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0.001 0.333 0.998 0.333],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0.001 0.001 0.998 0.333],'xtick',[],'ytick',[],'box','on','handlevisibility','off')


h2 = axes('Position',[0.2 0.333+0.1 0.75 0.20],'box','off');
par = 1;
x = power_vec(startfile:end);
yprime = allentries(startfile:end,par+1);
errorbar(x,yprime,allentries(startfile:end,par+10),'-bo','LineWidth',2);
ylabel('n (1e13/cm^2)')
xlabel('Power Density (kW/cm^2)')
xlim([2 5.75]) 
set(gca,'FontSize',fontsize,'box','off');

h1 = axes('Position',[0.2 0.667+0.1 0.75 0.20],'box','off');
par = 2;
errorbar(power_vec(startfile:end),allentries(startfile:end,par+1),allentries(startfile:end,par+10),'-bo','LineWidth',2);
ylabel('Temperature (K)')
xlabel('Power Density (kW/cm^2)')
xlim([2 5.75]) 
set(gca,'FontSize',fontsize,'box','off');

h3 = axes('Position',[0.2 0.1 0.75 0.20],'box','off');
par = 3;
hold on
errorbar(power_vec(startfile:end),allentries(startfile:end,par+1),allentries(startfile:end,par+10),'-bo','LineWidth',2);
ylabel('BGR (eV)')
xlabel('Power Density (kW/cm^2)')
xlim([2 5.75]) 
set(gca,'FontSize',fontsize,'box','off');

dim = [0.007 0.667+0.23 0.1 0.1];
str = {'(a)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.007 0.333+0.23 0.1 0.1];
str = {'(b)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.007 0.000+0.23 0.1 0.1];
str = {'(c)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
hold off

saveas(gcf,'fig_pars.png')

%%

%t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
pad = 0.06;
offset = 0.03;

fig4 = figure('Position',[0 0 1100 800]);

hold on

axes('Position',[0 0 0.5 0.5],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0 0.5 0.5 0.5],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0.5 0 0.5 0.5],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0.5 0.5 0.5 0.5],'xtick',[],'ytick',[],'box','on','handlevisibility','off')

h1 = axes('Position',[pad+offset 0.5+pad+offset 0.5-2*pad 0.5-2*pad]);
par = 4;
errorbar(power_vec(startfile:end),allentries(startfile:end,par+1),allentries(startfile:end,par+10),'-bo','LineWidth',2);
ylabel('Alpha (arb.)')
xlabel('Power Density (kW/cm^2)')
xlim([2 5.75]) 
set(gca,'FontSize',fontsize,'box','off');

h2 = axes('Position',[0.5+pad+offset 0.5+pad+offset 0.5-2*pad 0.5-2*pad]);
par = 5;
errorbar(power_vec(startfile:end),allentries(startfile:end,par+1),allentries(startfile:end,par+10),'-bo','LineWidth',2);
ylabel('Beta (arb.)')
xlabel('Power Density (kW/cm^2)')
xlim([2 5.75]) 
set(gca,'FontSize',fontsize,'box','off');

h3 = axes('Position',[pad+offset pad+offset 0.5-2*pad 0.5-2*pad]);
par = 7;
errorbar(power_vec(startfile:end),allentries(startfile:end,par+1),allentries(startfile:end,par+10),'-bo','LineWidth',2);
ylabel('C (arb.)')
xlabel('Power Density (kW/cm^2)')
xlim([2 5.75]) 
set(gca,'FontSize',fontsize,'box','off');

h4 = axes('Position',[0.5+pad+offset pad+offset 0.5-2*pad 0.5-2*pad]);
par = 8;
errorbar(power_vec(startfile:end),allentries(startfile:end,par+1),allentries(startfile:end,par+10),'-bo','LineWidth',2);
ylabel('D (arb.)')
xlabel('Power Density (kW/cm^2)')
xlim([2 5.75]) 
set(gca,'FontSize',fontsize,'box','off');
dim = [0.007 0.89 0.1 0.1];
str = {'(a)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.507 0.89 0.1 0.1];
str = {'(b)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.007 0.39 0.1 0.1];
str = {'(c)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.507 0.39 0.1 0.1];
str = {'(d)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');

saveas(gcf,'fig_alphabetaCD.png')




%% Fermi Levels, E_Thermal


Z = zeros(21,5);

for i = startfile:21
    
    n = allentries(i,2);
    T = allentries(i,3);
    
    if T > 294
        strain = polyval(f1,T);
    else
        strain = polyval(f1,294);
    end
    E_off_h = polyval(f3,strain); %eV
    m_e_K = polyval(f5,strain)*m0; %effective masses (highest fluence)
    m_h_K = polyval(f4,strain)*m0;
    m_h_G = polyval(f2,strain)*m0; 
    g_e_K = (degen_e_K*m_e_K)/(pi*hbar^2); %density of states
    g_h_K = (degen_h_K*m_h_K)/(pi*hbar^2);
    g_h_G = (degen_h_G*m_h_G)/(pi*hbar^2);

    const = [E_opt, E_bind, k_b, hbar, m0];

    E_fermi_e = F1(T,n.*1e13,const,g_e_K);
    E_fermi_h = F2(T,n.*1e13,const,g_h_K,g_h_G,E_off_h);

    a = 1e-4; %units

    n_h_K = a*g_h_K*k_b*T*log(exp(E_fermi_h/(k_b*T)) + 1);

    n_h_G = a*g_h_G*k_b*T*log(exp((E_fermi_h + E_off_h)/(k_b*T)) + 1);

    KKgap = @(T) (-353e-6)*T + 1.9934; %2.0118; %from Avinash calculations (E_thermal)
    E_thermal = E_opt - KKgap(T);

    Z(i,:) = [E_fermi_e, E_fermi_h, n_h_K*1e-13, n_h_G*1e-13, E_thermal];
end

%% Plots

figure();
set(gca,'FontSize',fontsize)
hold on
plot(power_vec(startfile:end),Z(startfile:end,1),'-ro','LineWidth',2);
plot(power_vec(startfile:end),Z(startfile:end,2),'-bo','LineWidth',2);
legend('Electron Ef','Hole Ef','Location','East')
xlabel('Power Density (kW/cm^2)')
ylabel('Energy (eV)')
hold off
saveas(gcf,'fermilevelvfluence.png')

%%

par = 6;
QE = allentries(startfile:end,par+1);
QE_errbar = allentries(startfile:end,par+10);
scale = 1/max(QE);
QE = scale*QE;
QE_errbar = scale*QE_errbar;


figure('Renderer', 'painters', 'Position', [10 10 650 800])
axes('Position',[0.001 0.001 0.998 0.498],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
axes('Position',[0.001 0.5 0.998 0.498],'xtick',[],'ytick',[],'box','on','handlevisibility','off')

h1 = axes('Position',[0.13 0.6 0.8 0.34]);
set(gca,'FontSize',fontsize)
hold on
plot(power_vec(startfile:end),Z(startfile:end,3),'-ro','LineWidth',2);
plot(power_vec(startfile:end),Z(startfile:end,4),'-bo','LineWidth',2);
plot(power_vec(startfile:end),allentries(startfile:end,2),'-ko','LineWidth',2);
xlabel('Power Density (kW/cm^2)')
ylabel('n_h (1e13/cm^2)')
legend('K-Valley','\Gamma Valley', 'Total','Location','East')
hold off


h2 = axes('Position',[0.13 0.1 0.8 0.34]);
hold on
%errorbar(power_vec(startfile:end),allentries(startfile:end,par+1),allentries(startfile:end,par+10),'-o','LineWidth',2);
errorbar(power_vec(startfile:end),QE,QE_errbar,'-bo','LineWidth',2);
set(gca,'FontSize',fontsize)
%plot(power_vec(startfile:end),peakheight(startfile:21),'r^','LineWidth',2);
ylabel('|\mu|^2 (arb.)')
xlabel('Power Density (kW/cm^2)')
xlim([2 5.75]) 
set(gca,'FontSize',fontsize)
dim = [0.007 0.895 0.1 0.1];
str = {'(a)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
dim = [0.007 0.395 0.1 0.1];
str = {'(b)'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20,'EdgeColor','none');
hold off


saveas(gcf,'fig_npervalley+dipole.png')

%%

alpha = allentries(:,5);
C = allentries(:,8);
beta = allentries(:,6);
D = allentries(:,9);


gam_e = alpha.*(Z(:,1).^2) + C;
gam_h = beta.*(Z(:,2).^2) + D;

tau_e = hbar./gam_e;
tau_h = hbar./gam_h;

figure();
hold on
plot(power_vec(startfile:end),tau_e(startfile:end))
plot(power_vec(startfile:end),tau_h(startfile:end))
ylim([0 1e-12])
hold off

%%

% Y = allscans2D;
% 
% % allscans2D = allscans2D(allscans2D(:,3) > 301,:);
% 
% Y = zeros(files,8);
% for i = 1:files
%     scan_list = allscans2D(allscans2D(:,1) == i,:);
%     scan_list = sortrows(scan_list,8);
%     try
%         Y(i,:) = scan_list(1,:);
%     catch
%     end
% end


%% BGR check

endfile = 15;

x = allentries(startfile:endfile,2);
z = allentries(startfile:endfile,4);

a_bohr = 6.0e-8;
bgr = @(n) 3.1*(((n*1e13).*a_bohr^2).^(1/3)).*(E_bind);
nspace = linspace(2.5,5,1000);
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
plot(z,x,'*','MarkerSize',8)
plot(bgrvals,nspace,'-r','LineWidth',2);
ylim([2 5.5])
xlabel('Bandgap Renormalization Energy (eV)')
ylabel('Charge Carrier Density (1e13/cm^2)')
legend('Fit Results','Predicted BGR vs Charge Carrier Density','Location','SouthEast')
set(gca,'FontSize',fontsize)
hold off

saveas(gcf,'fig_nvsBGR.png')


%%

xprime = logspace(-3,1);

yprime = -3.1*xprime.^(1/3);
y2prime = -2.8*xprime.^(1/3);
y3prime = -2.5*xprime.^(1/3);

E0 = E_bind; %30.1588;

na_0squared = x*10^13*a_bohr^2;
BGR_norm = -z./E0;

semilogx(xprime,yprime,'-',xprime,y2prime,'-',xprime,y3prime,'-',na_0squared,BGR_norm,'*')
legend('3.1','2.8','2.5','data')
ylabel('E_{BGR}/E_{bind}')
xlabel('na_{b}^2')
ylim([-2 0])
set(gca,'FontSize',fontsize)




