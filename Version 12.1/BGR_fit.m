

%p0 = [3.1, 0.42];
p0 = [3.1, 6.5e-8, 0.42];

g = fittype( @(a,b,c,x) a*(((x*1e13).*b^2).^(1/3)).*c );

x = allentries(startfile:end,2);
z = allentries(startfile:end,4);

curve = fit( x, z, g, 'StartPoint', p0  ,'Robust', 'LAR' );
disp(curve)

hold on
plot(curve,x,z,'o')
legend('Location','NorthWest')
hold off