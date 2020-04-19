function Ef_h = F2(T,n,const,g_h_K,g_h_G,E_off_h)

    a = 1e-4; %units

    k_b = const(3);

    p_fun = @(T,Ef_h) a*k_b*T*g_h_K*log(exp(Ef_h/(k_b*T)) + 1) + ...
        a*k_b*T*g_h_G*log(exp((Ef_h + E_off_h)/(k_b*T)) + 1);
    
    points = 1000; %resolution of interpolation

    ef_list = linspace(-0.5,0.5,points);

    v2 = p_fun(T,ef_list);
    
    Ef_h = interp1(v2,ef_list,n);
end