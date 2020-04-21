
clear all
close all

filenum = 14;


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
tic

global ticker finish filenum

ticker = 6;
finish = 0;

while finish == 0
    

        disp('Starting')
        list = fliplr(1:21);
        filenum = list(ticker);
        
        filestring = num2str(filenum);
        filename = strcat('montecarloerrorsaveVersion12.11_April1720_file',filestring,'.mat');
        S2 = load(filename);
        
        disp('Loaded File')
        
        X = S2.X;
        disp(X)
        disp(S2.RESNORM)
        %LSQ_FIT_MAIN_SCRIPT %get fit X
        
        if S2.RESNORM < 0.10
            
            disp('Valid X')
            try
                Monte_Carlo_Error %get X_err 250 trials
                save(filename)
            catch
                disp('Error')
            end
            ticker = ticker+1;
            disp('Moving to Next File...')
            
        end
           
        clearvars -except ticker finish
        close all
        
        if ticker > 21
            finish = 1;
        end

end

toc