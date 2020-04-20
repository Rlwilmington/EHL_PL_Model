function out = D_e(const,E,n,T,alpha,C,f1,f5)

degen_e_K = 2; %band degeneracy

k_b = const(3);
hbar = const(4);
m0 = const(5);

if T > 294
    strain = polyval(f1,T);
else
    strain = polyval(f1,294);
end

m_e_K = polyval(f5,strain)*m0; %effective masses (highest fluence)
g_e_K = (degen_e_K*m_e_K)/(pi*hbar^2); %density of states

Ef_e = F1(T,n,const,g_e_K);

conv_e = @(tau) g_e_K.*( (exp((tau - Ef_e)./(k_b.*T)) + 1).^-1 ) .*...
    (1./(2*pi)).*((C + alpha*(tau-Ef_e).^2)./((E - tau).^2 + ((C + alpha*(tau-Ef_e).^2)./2).^2));


out = integral(conv_e,0,Inf,'ArrayValued',true); %zero to inf means g_e is only evaluated for positive values, so the density of states is appropriately valued for only energy values outside the gap

end

