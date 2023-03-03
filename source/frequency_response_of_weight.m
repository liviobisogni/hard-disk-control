%
% Il seguente script plotta le risposte in frequenza del Modello del
% sistema e delle varie funzioni peso utili alla sintesi del controllore.

% omega = logspace(-2,2,200);

close all; clc;

global save_figure_flag;

% Carico modello, funzioni peso, ecc
olp_hdd;

wts_hdd;


% frequency response of G       %%SUPERFLUO??
olp_ic = G(2);
w = logspace(2,5,500);
figure(1)
olp_ic30 = usample(olp_ic,30);
bode(olp_ic.NominalValue,'r-',olp_ic30,'b--',w), grid
title('Bode plot of the uncertain plant')
legend('Nominal system','Random samples','Location','SouthWest') %legend('Nominal system','Random samples',3)
if (save_figure_flag == 1)
    fig_name = ['Bode_plot_of_the_uncertain_plant.png'];
    saveas(gcf, fig_name)
end

% Frequency responde of M (Wm)
w = logspace(1,7,500);
figure(2)
bodemag(Wm,w), grid
title('Frequency Responce of M (W_m, the Model Transfer Function)');
if (save_figure_flag == 1)
    fig_name = ['Frequency_Responce_of_M_(W_m).png'];
    saveas(gcf, fig_name)
end
% figure(1), 
% sigma(Wm(1,1), 'b-')
% grid
% xlabel('Frequency')
% ylabel('Magnitude')
% title('Frequency Responce of M (W_m, the Model Transfer Function)');

% Frequency responce of inverse of Wp
w = logspace(-1,5,500);
figure(3)
bodemag(1/Wp,w), grid
title('Inverse Performance Weighting Function (1/W_p)')
if (save_figure_flag == 1)
    fig_name = ['Inverse_Performance_Weighting_Function_(W_p)^-1.png'];
    saveas(gcf, fig_name)
end
% figure(2), 
% sigma(inv(Wp(1,1)),'b-')
% grid
% xlabel('Frequency')
% ylabel('Magnitude')
% title('Frequency Responce of the inverse of W_p (the Inverse Performance Weighting Function)');


%Frequency responce of Wn
w = logspace(0,4,500);
figure(4)
bodemag(Wn,'r',w), grid
title('Noise shaping filter frequency response (W_n)')
if (save_figure_flag == 1)
    fig_name = ['Noise_shaping_filter_frequency_response_(W_n).png'];
    saveas(gcf, fig_name)
end
% figure(3), 
% sigma(Wn(1,1),'b-')
% grid
% xlabel('Frequency')
% ylabel('Magnitude')
% title('Frequency Responce of W_n (the Noise Shaping Function)');

%Frequency responce of Wu
w = logspace(-1,6,500);
figure(5)
bodemag(Wu,'r',w), grid
title('Control Weighting Function (W_u)')
if (save_figure_flag == 1)
    fig_name = ['Control_Weighting_Function_(W_u).png'];
    saveas(gcf, fig_name)
end
% figure(4), 
% sigma(Wu(1,1),'b-')
% grid
% xlabel('Frequency')
% ylabel('Magnitude')
% title('Frequency Responce of W_u (the Control Weighting Function)');
