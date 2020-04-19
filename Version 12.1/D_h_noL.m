function [out] = D_h_noL(const,E,n,T,f1,f2,f3,f4)

degen_h_K = 2;
degen_h_G = 1;

k_b = const(3);
hbar = const(4);
m0 = const(5);

if T > 293
    strain = polyval(f1,T);
else
    strain = 0.6934;
end

E_off_h = polyval(f3,strain); %eV
m_h_K = polyval(f4,strain)*m0;
m_h_G = (f2(1)*exp(-f2(2)*strain) + f2(3))*m0;
g_h_K = (degen_h_K*m_h_K)/(pi*hbar^2);
g_h_G = (degen_h_G*m_h_G)/(pi*hbar^2);

Ef_h = F2(T,n,const,g_h_K,g_h_G,E_off_h);

% conv_h = @(tau) ( g_h_K.*( (exp((E - tau - Ef_h)./(k_b.*T)) + 1).^-1 ) + g_h_G.*( (exp((E - tau - Ef_h - E_off_h)./(k_b.*T)) + 1).^-1 ));
% 
% out = integral(conv_h,0,E,'ArrayValued',true);
%

out = (g_h_K.*( (exp((E - Ef_h)./(k_b.*T)) + 1).^-1 ) + g_h_G.*( (exp((E - Ef_h - E_off_h)./(k_b.*T)) + 1).^-1 ));

end

