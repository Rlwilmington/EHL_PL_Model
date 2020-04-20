
clear all
close all

filenum = 21;


LSQ_FIT_MAIN_SCRIPT %get fit X

%Monte_Carlo_Error %get X_err 250 trials

filestring = num2str(filenum);
filename = strcat('montecarloerrorsaveVersion12.11_April1720_file',filestring,'.mat');

save(filename)




%%

% for i = 1:21
%     
% 
% 
%         %list = fliplr(1:21);
%         filenum = i;
% 
%         LSQ_FIT_MAIN_SCRIPT %get fit X
% 
%         %Monte_Carlo_Error %get X_err 250 trials
% 
%         filestring = num2str(filenum);
%         filename = strcat('montecarloerrorsaveVersion12.10_April1720_file',filestring,'.mat');
% 
%         save(filename)
%         
%         clear all
%         close all
% 
% end

%%


global i finish

i = 1;
finish = 0;

while finish == 0
    


        %list = fliplr(1:21);
        filenum = i;
        
        filestring = num2str(filenum);
        filename = strcat('montecarloerrorsaveVersion12.11_April1720_file',filestring,'.mat');
        S = load(filename);
        
        %LSQ_FIT_MAIN_SCRIPT %get fit X
        
        if S.RESNORM < 0.10
        

            Monte_Carlo_Error %get X_err 250 trials
            save(filename)
        
            i = i+1;
            
        end
           
        clearvars -except i finish
        close all
        
        if i > 21
            finish = 1;
        end

end
