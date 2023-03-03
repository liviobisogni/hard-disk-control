% H_infinity design for the Disk Drive System
%

close all; clc;

global K; global save_figure_flag;

% Carico modello, funzioni peso, ecc
olp_hdd;

% Carico i flag
weight_frsp_flag = evalin('base','weight_frsp_flag');

%%
% Visualizzo risposte in frequenza delle funzioni peso
if ( weight_frsp_flag == 1 )
        frequency_response_of_weight;
end


nmeas = 1;
ncon = 1;
gmin = 0;
gmax = 10;
tol = 0.001;
hin_ic = sys_ic.Nom;
[K_hin,clp_hin,gopt] = hinfsyn(hin_ic,nmeas,ncon,gmin,gmax,tol);
gmin = 1.1*gopt;
[K_hin,clp_hin,gfin] = hinfsyn(hin_ic,nmeas,ncon,gmin,gmin,tol);
disp(' ')
get(K_hin)
disp(' ')
disp('Closed-loop poles')
%sp = pole(clp_hin) %CANC?
%clp_poles = pole(clp_hin)
omega = logspace(-2,6,100);
figure
sigma(clp_hin,omega)
title('Singular Value Plot of clp')
xlabel('Frequency (rad/sec)')
ylabel('Magnitude')
if (save_figure_flag == 1)
    fig_name = ['K_hinf__Singular_Value_Plot_of_clp.png'];
    saveas(gcf, fig_name)
end
K = K_hin;
clp_ic = lft(sys_ic,K,1,1);
clp_poles = pole(clp_ic)


figure,
%sigma(K_hin,logspace(1,5,30)), grid
bodemag(K_hin,logspace(1,7,500)), grid
%xlabel('Frequency (rad/min)'), ylabel('Magnitude')
title('Frequency response of the controller')
if (save_figure_flag == 1)
    fig_name = ['K_hinf__Frequency_response_of_the_controller.png'];
    saveas(gcf, fig_name)
end