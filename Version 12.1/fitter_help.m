
x = linspace(0,2.5,100);

f = @(x) 3*exp(-0.5*x) + 0;

hold on
plot(x,f(x),'x')
plot(masses.Strain,masses.G_VB,'*')
hold off