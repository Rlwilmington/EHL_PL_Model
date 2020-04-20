function [out] = D_h(const,E,n,T,beta,D,f1,f2,f3,f4)

degen_h_K = 2;
degen_h_G = 1;

k_b = const(3);
hbar = const(4);
m0 = const(5);

if T > 294
    strain = polyval(f1,T);
else
    strain = polyval(f1,294);
end

E_off_h = polyval(f3,strain); %eV
m_h_K = polyval(f4,strain)*m0;
m_h_G = polyval(f2,strain)*m0;
g_h_K = (degen_h_K*m_h_K)/(pi*hbar^2);
g_h_G = (degen_h_G*m_h_G)/(pi*hbar^2);

Ef_h = F2(T,n,const,g_h_K,g_h_G,E_off_h);

conv_h = @(tau) ( g_h_K.*( (exp((tau - Ef_h)./(k_b.*T)) + 1).^-1 ) + g_h_G.*( (exp((tau - Ef_h - E_off_h)./(k_b.*T)) + 1).^-1 ))...
    .*(1./(2*pi)).*((D + beta*(tau-Ef_h).^2)./((E - tau).^2 + ((D + beta*(tau-Ef_h).^2)./2).^2));

out = integral(conv_h,0,Inf,'ArrayValued',true);

end

