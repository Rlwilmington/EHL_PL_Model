function [X, RESNORM] = PL_fit_generator(par,const,xdata,ydata,f1,f2,f3,f4,f5,maxevals)

    minimizer = @(par) PL_resid(par,const,xdata,ydata,f1,f2,f3,f4,f5);

    %par0 = [A,n,T,E_BGR,alpha,beta]; %init

    UB = [6,  700,  0.8,  100, 100,   1, 0.05, 0.05]; %define fit bounds
    LB = [0.3,  280,  0.4,  0,  0, 0, -0.01, -0.01];

    fitops = optimoptions('lsqnonlin','maxfunevals',maxevals,'TolFun',1e-6,'TolX',1e-6,'Algorithm','trust-region-reflective');

    %fitops = optimoptions('lsqnonlin','maxfunevals',maxevals,'TolFun',1e-9,'TolX',1e-9,'Display','iter-detailed','Algorithm','trust-region-reflective');
    
    %options.Algorithm = 'levenberg-marquardt'; %select other algorithm

    [X, RESNORM] = lsqnonlin(minimizer, par, LB, UB, fitops);

end


