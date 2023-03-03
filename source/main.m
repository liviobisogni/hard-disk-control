%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main file of the Robust Control of a Hard Disk Drive Poject %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Definisco variabili globali
% Init_global;

close all; clear all; clc; closereq;  % pulisco
%mod_hdd     % crea il modello del sistema incerto
olp_hdd;    % crea il modello del sistema a ciclo aperto incerto

% Controllori (definiti come variabili globali)
global K; global K_hin; global K_mu;global K_lsh; global K_PID; global K_lqgltr_u; global K_lqgltr_y; global CL_Nominal;
global save_figure_flag;
save_figure_flag = 0;       % porre uguale a 1 solo se si vuole salvare le figure

% Flag = 1 --> Visualizza le risposte in frequenza delle funzioni peso
weight_frsp_flag = 0;


GUI;    % apre l'interfaccia grafica