function res_y = PL_resid(par,const,xdata,ydata,f1,f2,f3,f4,f5)

    [~, xfit_conv, yfit, ~, ~] = PLfitfun_int(par,const,f1,f2,f3,f4,f5);
    
    yfit_interp = interp1(xfit_conv,yfit,xdata);
    
    res_y = yfit_interp - ydata;

end