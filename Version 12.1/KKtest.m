KKgap1 = @(T) (-353e-6)*(T)+1.993429; %from Avinash calculations (E_thermal)
KKgap2 = @(T) (353e-6)*(T-293); %from Avinash calculations (E_thermal)

Temps = linspace(180,650,100);
E_T1 = zeros(100,1);
E_T2 = zeros(100,1);


for i = 1:100
    
    T = Temps(i);
    
    if T > 293
        E_thermal1 = E_opt - KKgap1(T);
        E_thermal2 = KKgap2(T);
    else
        E_thermal1 = 0;
        E_thermal2 = 0;
    end

    E_T1(i) = E_thermal1;
    E_T2(i) = E_thermal2;
    
end
