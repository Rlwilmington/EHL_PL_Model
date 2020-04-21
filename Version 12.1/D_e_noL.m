function out = D_e_noL(const,E,n,T,f1,f5)

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

%conv_e = @(tau) g_e_K.*( (exp((E - tau - Ef_e)./(k_b.*T)) + 1).^-1 );


%out = integral(conv_e,0,E,'ArrayValued',true); %zero to inf means g_e is only evaluated for positive values, so the density of states is appropriately valued for only energy values outside the gap

out = g_e_K.*( (exp((E - Ef_e)./(k_b.*T)) + 1).^-1 );

end

