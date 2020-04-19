function [f1, f2, f3, f4, f5] = StrainvTemp_MAC()

%% Load Data


%Mac OS X Version
load LatticeExpVSPower.mat LatticeExpVSPower 
load TempVsPower.mat TempVsPower




%%

% M1 = sortrows(LatticeExpVSPower,1);
% M2 = sortrows(TempVsPower,1);
% 
% fit1 = fit(M2(:,1),M1(:,1),'poly1');
% f1 = [fit1.p1, fit1.p2];


%% Calibrate Sample Temperature to Percent Strain Using Raman Data

LatticeExpVSPower(:,1) = flipud(LatticeExpVSPower(:,1));
tempStrainCoeff = polyfit(TempVsPower(1:6,1),LatticeExpVSPower(1:6,1),1);
f1 = tempStrainCoeff;

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

f3 = polyfit(M(:,1),M(:,4),1);

%%

masses = readtable('MoS2_masses.txt','Delimiter',' ','ReadVariableNames',false);
masses.Properties.VariableNames = {'Strain','LatExp','K_CBL','K_CBU','K_VBU','K_VBL','G_VB','Q_CBL','Q_CBU'};

% massInterp = table;
massInterp = table('Size',[numPts 5],'VariableTypes',{'double','double','double','double','double'},'VariableNames',{'Strain','K_CB','K_VB','Q_CB','G_VB'});

massInterp.Strain = latticeExp.';
massInterp.K_CB = interp1(masses.Strain,masses.K_CBL,latticeExp,'linear')';
massInterp.K_VB = interp1(masses.Strain,masses.K_VBU,latticeExp,'linear')';
massInterp.G_VB = interp1(masses.Strain,masses.G_VB,latticeExp,'linear')';

%g = fittype('a*exp(-b*x) + c');

%fit2 = fit(masses.Strain,masses.G_VB,g,'Start',[3, 0.5, 0]); %effective mass G-band holes (m_h_G)
%f2 = [fit2.a, fit2.b, fit2.c];
f2 = polyfit(masses.Strain,masses.G_VB,1);

f4 = polyfit(masses.Strain,masses.K_VBL,1); %effective mass K-band holes (m_h_K)

f5 = polyfit(masses.Strain,masses.K_CBL,1); %effective mass K-band electrons (m_e_K)


%% Plots
% 
% x1 = linspace(300,600);
% y1 = polyval(f1,x1);
% 
% figure();
% hold on
% plot(TempVsPower(1:6,1),LatticeExpVSPower(1:6,1),'*')
% plot(x1,y1,'-')
% hold off
% 
% x2 = linspace(0,3,100);
% y2 = polyval(f2,x2);
% 
% figure();
% hold on
% plot(masses.Strain,masses.G_VB,'*')
% plot(x2,y2,'-')
% hold off
% 
% x3 = linspace(0,2.5,100);
% y3 = polyval(f3,x3);
% 
% figure();
% hold on
% plot(M(:,1),M(:,4),'*')
% plot(x3,y3,'-')
% hold off
% 
% x4 = linspace(0,3,100);
% y4 = polyval(f4,x4);
% 
% figure();
% hold on
% plot(masses.Strain,masses.K_VBL,'*')
% plot(x4,y4,'-')
% hold off
% 
% x5 = linspace(0,3,100);
% y5 = polyval(f5,x5);
% 
% figure();
% hold on
% plot(masses.Strain,masses.K_CBL,'*')
% plot(x5,y5,'-')
% hold off
% 
% 


end

