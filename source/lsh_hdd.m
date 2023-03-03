% Loop Shaping Design for the Hard Disk Drive System
%

close all; clc;

global K; global save_figure_flag;

% Carico modello, funzioni peso, ecc
olp_hdd;

mod_hdd
%
% set the precompensator
nuW1 = [0.05 1];
dnW1 = [1 0];
gainW1 = 4;
W1 = gainW1*tf(nuW1,dnW1);
%
% compute the loop shaping controller
[K_0,cl,gam,info] = ncfsyn(G(1,2).Nom,W1);
%clp_ploes = pole(cl)
emax = info.emax;
disp(['The nugap robustness emax = ' num2str(emax)]);
%
% frequency responses of the plant and shaped plant
Gs = info.Gs;
w = logspace(-1,5,400);
figure(1)
bodemag(G(1,2).Nom,'b-',Gs,'r--',w), grid
title('Frequency responses of the plant and shaped plant')
legend('Original plant','Shaped plant')
if (save_figure_flag == 1)
    fig_name = ['K_lsh__Frequency_responses_of_the_plant_and_shaped_plant.png'];
    saveas(gcf, fig_name)
end

% obtain the negative feedback controller
K_lsh = -K_0;
K = K_lsh;
clp_ic = lft(sys_ic,K,1,1);
clp_poles = pole(clp_ic)

figure,
%sigma(K_lsh, logspace(-3,3,100)), grid
bodemag(K_lsh,logspace(1,7,500)), grid
%xlabel('Frequency (rad/min)'), ylabel('Magnitude')
title('Frequency response of the controller')
if (save_figure_flag == 1)
    fig_name = ['K_lsh__Frequency_response_of_the_controller.png'];
    saveas(gcf, fig_name)
end