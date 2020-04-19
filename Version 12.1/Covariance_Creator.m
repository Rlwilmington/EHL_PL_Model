

maxevals = 2000;
t = readtable('saverunFeb1620.txt');
allentries = table2array(t);

N = length(allentries);
xdata = EHL(:,1);

UB = [6,  700,  0.8,  100, 100,   1, 0.05, 0.05]; %define fit bounds
LB = [0.3,  280,  0.4,  0,  0, 0, -0.01, -0.01];

fitops = optimoptions('lsqnonlin','maxfunevals',maxevals,'TolFun',1e-9,'TolX',1e-9,...
    'disp','iter-detailed','Algorithm','trust-region-reflective');

dof = length(xdata) - length(X);

for i = 2:N
    
    filenum = allentries(i,1);
    
    ydata = 0.1*data.y(:,filenum);
    
    par = allentries(i,2:9);
    
    minimizer = @(par) PL_resid(par,const,xdata,ydata,f1,f2,f3,f4,f5); %define single-var residual function to be minimized

    [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(minimizer, par, LB, UB, fitops); %perform LSQ fit
   
    J = JACOBIAN;
    J = full(J);
    mse = (1/dof)*sum(RESIDUAL.^2); %variance of residuals (assuming ~uniform)
    cov = (J.'*J)*mse;
    SE = sqrt(diag(cov)).';
        
    entry = horzcat(filenum,X,RESNORM,SE);

    entry_log = vertcat(entry_log,entry);
    
    disp(entry(1:10))
    
end

Tlog = array2table(entry_log);
writetable(Tlog,'entrylogVersion10withCovFeb2820.txt');


