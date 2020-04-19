
clear all
close all

filenum = 21;


LSQ_FIT_MAIN_SCRIPT %get fit X

%Monte_Carlo_Error %get X_err 250 trials

filestring = num2str(filenum);
filename = strcat('montecarloerrorsaveVersion12.8_April1720_file',filestring,'.mat');

save(filename)




%%

for i = 1:21
    


        list = fliplr(1:21);
        filenum = list(i);

        LSQ_FIT_MAIN_SCRIPT %get fit X

        %Monte_Carlo_Error %get X_err 250 trials

        filestring = num2str(filenum);
        filename = strcat('montecarloerrorsaveVersion12.8_April1720_file',filestring,'.mat');

        save(filename)
        
        clear all
        close all

end


