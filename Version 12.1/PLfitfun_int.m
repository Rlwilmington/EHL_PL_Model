function [xfit_g, xfit_conv, yfit, u_vec, v_vec] = PLfitfun_int(par,const,f1,f2,f3,f4,f5)


    E_opt = const(1);
    E_bind = const(2);

    %disp(par)
    %2.0118; %from Avinash calculations (E_thermal)%from Avinash calculations (E_thermal)
    T = par(2);
    
    KKgap = @(T) (353e-6)*(T-294);
    
    if T > 294
        E_thermal = KKgap(T);
    else
        E_thermal = 0;
    end
    
    n = 1e13*par(1); %already converted to 1/m^2
    T = par(2);
    E_BGR = par(3);
    alpha = par(4); %Lorentzian 'gamma' -> alpha*(E-Ef)^2 + beta
    beta = par(5);
    A = 1e-38*par(6);
    C = par(7);
    D = par(8);
    
    %%%% CREATE NEW X REGION %%%%
    step = 0.005;
    xfit_g = (-1:step:1).'; %end at 4eV, though this choice should not make a difference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %% Add in Effmass, Eoff dep on Temp
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     strain = f1(T);
%     E_off_h = f3(strain); %eV
%     m_e_K = f5(strain)*m0; %effective masses (highest fluence)
%     m_h_K = f4(strain)*m0;
%     m_h_G = f2(strain)*m0; %1.9*m0;
%     g_e_K = (degen_e_K*m_e_K)/(pi*hbar^2); %density of states
%     g_h_K = (degen_h_K*m_h_K)/(pi*hbar^2);
%     g_h_G = (degen_h_G*m_h_G)/(pi*hbar^2);
%     
%     const = [E_opt, E_bind, k_b, g_e_K, g_h_K, g_h_G, hbar, E_off_h, E_thermal, m0]; %constants necessary for fit model construction
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 
    u = @(E) D_e(const,E,n,T,alpha,C,f1,f5);
    v = @(E) D_h(const,E,n,T,beta,D,f1,f2,f3,f4);

    yfit = A*conv(u(xfit_g),v(xfit_g),'same'); %perform numerical convolution
    u_vec = (1/max(u(xfit_g)))*u(xfit_g);
    v_vec = (1/max(v(xfit_g)))*v(xfit_g);
    
    xfit_conv = xfit_g + (E_opt + E_bind - E_thermal- E_BGR); %convert hnuprime to hnu
    
end