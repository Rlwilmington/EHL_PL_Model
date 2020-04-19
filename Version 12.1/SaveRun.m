

%% Bulk Save for Fits

t = readtable('SaveDataVersion9Feb282020.txt');
allentries = table2array(t);

entry = horzcat(filenum,X,RESNORM,SE.');

allentries = vertcat(allentries,entry);

allentries = sortrows(allentries,1);
T = array2table(allentries);
writetable(T,'SaveDataVersion9Feb282020.txt')

disp(allentries(:,1:10))


%% Meta Outputs

Yprime = sortrows(allentries_2,1);
Y = Yprime(2:end,:); %cut off init row

figure();
semilogy(Y(:,1),Y(:,8),'*')
title('Fit Quality')
ylabel('MSE')
xlabel('File #')

a_bohr = 6.5e-8;
bgr = @(n) 3.1*(((n*1e13).*a_bohr^2).^(1/3)).*(E_bind);
nspace = linspace(0.4,3,1000);
bgrvals = bgr(nspace);

figure();
hold on
plot(nspace,bgrvals,'r')
plot(Y(:,2),Y(:,4),'*')
hold off
title('Number Density vs Bandgap Renormalization Energy')
xlabel('Number Density (10^13 carriers/cm^2)')
ylabel('BGR (eV)')

figure();
plot(Y(:,1),Y(:,2),'*')
title('Number Density vs Fluence')

figure();
plot(Y(:,1),Y(:,3),'*')
title('Temperature vs Fluence')

figure();
plot(Y(:,1),Y(:,4),'*')
title('BGR vs Fluence')

figure();
plot(Y(:,1),Y(:,5),'*')
title('Alpha vs Fluence')

figure();
plot(Y(:,1),Y(:,6),'*')
title('Beta vs Fluence')

figure();
semilogy(Y(:,3),Y(:,8),'*')
title('Temperature vs Fit Quality')













