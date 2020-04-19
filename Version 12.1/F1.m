function Ef_e = F1(T,n,const,g_e_K)

    a = 1e-4; %units

    k_b = const(3);

    n_fun = @(T,Ef_e) a*k_b*T*g_e_K*log(exp(Ef_e/(k_b*T)) + 1);
    
    points = 1000; %resolution of interpolation

    ef_list = linspace(-0.5,0.5,points);

    v1 = n_fun(T,ef_list);
    
    Ef_e = interp1(v1,ef_list,n);
end