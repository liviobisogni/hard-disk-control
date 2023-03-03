% mu-synthesis of the Disk Drive System
%


close all; clc

global K; global save_figure_flag;

% Carico modello, funzioni peso, ecc
olp_hdd;

% Carico i flag
%weight_frsp_flag = evalin('base','weight_frsp_flag');

%%
% % Visualizzo risposte in frequenza delle funzioni peso
% if ( weight_frsp_flag == 1 )
%         frequency_response_of_weight;
% end

nmeas = 1;
ncont = 1;
mu_ic = usubs(sys_ic,'w1','nom','w2','nom','w3','nom','w4','nom', ...
                     'z1','nom','z2','nom','z3','nom','z4','nom')
fv = logspace(1,5,30);
opt = dkitopt('FrequencyVector',fv, ...
               'DisplayWhileAutoIter','on', ...
              'NumberOfAutoIterations',4)
[K_mu,CL_mu,bnd_mu,dkinfo] = dksyn(mu_ic,nmeas,ncont,opt);
K = K_mu;
%clp_poles = pole(CL_mu)
clp_ic = lft(sys_ic,K,1,1);
clp_poles = pole(clp_ic)



figure,
%sigma(K_mu,fv),grid
bodemag(K_mu,logspace(1,7,500)), grid
%xlabel('Frequency (rad/min)'), ylabel('Magnitude')
title('Frequency response of the controller')
if (save_figure_flag == 1)
    fig_name = ['K_mu__Frequency_response_of_the_controller.png'];
    saveas(gcf, fig_name)
end