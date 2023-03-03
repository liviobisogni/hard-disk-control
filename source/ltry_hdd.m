%
% Lo script provvede alla sintesi del controllore LQG/LTR con apertura del
% loop in ingresso al sistema

close all; clc; 

% Carico modello, funzioni peso, ecc
olp_hdd;

% % prefiltro G_shaped =  W1GW2
% nuW1 = [1.1   1];
% dnW1 = [10    0];
% gainW1 = 1.7;
% w1 = nd2sys(nuW1,dnW1,gainW1);
% W1 = daug(w1,w1);

% set the precompensator
nuW1 = [0.05 1];
dnW1 = [1 0];
gainW1 = 4;
W1 = gainW1*tf(nuW1,dnW1);


% compute the loop shaping controller
[K_0,cl,gam,info] = ncfsyn(G(1,2).Nom,W1);  %!!!
emax = info.emax;
%disp(['The nugap robustness emax = ' num2str(emax)]);

Gs = info.Gs;
%ss(Gs)
%get(Gs)
%disp('!!!!!!!')
G_shaped = Gs;
%ss(G_shaped)
%get(G_shaped)
% G_shaped = mmult(pck(A,B,C,D),W1);
[A,B,C,D] = unpck(G_shaped);
A;
B;
C;
D;
sys_shaped = ss(Gs.A,Gs.B,Gs.C,Gs.D);
%ss(sys_shaped)
%get(sys_shaped)

% Senza prefiltro
% basti commentare la parte sopra

%%
% Inizializzo matrici peso per LQR e KF
W = Gs.B*Gs.B'; % varianza eumore su stato
%W = B*B' % varianza eumore su stato
V = eye(1); % varianza rumore misura
%V = eye(2); % varianza rumore misura
Q = Gs.C'*Gs.C; % peso stati
%Q = C'*C % peso stati
R = eye(1); % peso controlli
%R = eye(2); % peso controlli
WV = blkdiag(W,V);
QR = blkdiag(Q,R);

% Frequenze per plot
omega = logspace(-5,5);

%% 
% Creo sostanzialmente LQG

% % Design LQR
% Kc = lqrc(A,B,QR);

% Design KF
Kf = (lqrc(Gs.A',Gs.C',WV))';
%Kf = (lqrc(A',C',WV))';

% Valori singolari di T_kf, ? il modello da recuperare T_kf = C(.)Kf
T_kf = ss(Gs.A,Kf,Gs.C,0);
%T_kf = ss(A,Kf,C,zeros(2,2));
figure,
sigma(T_kf, 'r-', omega), grid % Desired Loop Transfer Matrix Responce
xlabel('Frequency (rad/min)'), ylabel('Magnitude')
title('Desired Loop Transfer Matrix Responce')

% Valori di recupero
% q = [1e3 1e6 1e9 1e12];
% q = [0 1 1e2 1e3 1e6 1e9 1e12];

%% *********************************    Loop Transfer Recovery   ***********************************
%
%Sintesi mediante ltrsyn. Il controllore restituito ? quello per q(end)

% [K_lqgltr] = ltrsyn(sys_shaped, Kf, Q, R, q, omega, 'OUTPUT');
% hold on
% sigma(T_kf, 'r-.', omega); % Desired Loop Transfer Matrix Responce
% grid

% Recupero per 1e9
figure,
[K_lqgltr] = ltrsyn(sys_shaped, Kf, Q, R, q, omega, 'OUTPUT'); % q(3)
hold on
sigma(T_kf, 'r-.', omega),grid % Desired Loop Transfer Matrix Responce

%%
% % Sintesi manuale. Anche qui il controllore ? quello ottenuto per q(end)
% 
% figure(2)
% hold on
% sigma(T_kf, 'r-.', omega)
% st='T desired';
% 
% for (i = 1:length(q))
%     
%     % Incremento la matrice di peso degli stati
%     Q = Q + q(i)*C'*C; % Qc = Q_0 + q^2 (C'VC), Q_0 = H'*H
%     QR = blkdiag(Q,R);
%     % Calcolo il nuovo guadagno dell'LQR
%     Kc = lqrc(A,B,QR);
%     % Compongo il controllore
%     A_K = A - B*Kc - Kf*C + Kf*D*Kc;            B_K= -Kf;                                       %     A_K = A - B*Kc - Kf*C + Kf*D*Kc;          B_K= Kf;
%     C_K = Kc;                                   D_K = zeros(size(Kc)*[1;0],size(Kf)*[0;1]);
%     
%     K_lqgltr = pck(A_K,B_K,C_K,D_K);
%     
%     % E' la loop m.d.t che bisogna far tendere a T_kf
%     T_ltr = mmult(G_shaped, K_lqgltr); % E' la loop m.d.t che bisogna far tendere a T_lq
%     
%    % Plotto valori singolari di T_ltr 
%     sigma(T_ltr, omega)
%    
%     st =[st strcat(string('T LTR SV for q='),num2str(q(i)))];
%     legend(st(:))
%     
% end
% 
%  grid

%%
% % Compongo il controllore Klqgltr con il prefiltro 
% n_filter = [1.1   1];
% d_filter = [10    0];
% g_filter = 1.7;
% w_filter = g_filter*tf(n_filter,d_filter);
% W_filter = blkdiag(w_filter,w_filter);

% set the precompensator
n_filter = [0.05 1];
d_filter = [1 0];
g_filter = 4;
w_filter = g_filter*tf(n_filter,d_filter);
W_filter = w_filter;

K_lqgltr_y = W_filter*tf(K_lqgltr);
K_lqgltr_y = ss(K_lqgltr_y);
K=K_lqgltr_y;
clp_ic = lft(sys_ic,K,1,1);
clp_poles = pole(clp_ic)

figure,
sigma(K_lqgltr_y,logspace(-3,3,100)), grid
xlabel('Frequency (rad/min)'), ylabel('Magnitude')
title('Frequency responses of the controller')

clear n_filter; clear d_filter; g_filter; clear w_filter; clear W_filter;


