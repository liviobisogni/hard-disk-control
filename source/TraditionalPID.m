%mio

%function [K,Tcl] = TraditionalPID( G, WC, OS, AMP, ST )
%function [K,Tcl] = TraditionalPID( G, WC, OS, AMP, ST, GM, PM )

%>> TraditionalPID(G,628e3,20,20,1.5)
%>> showTunable(CL_final)

clc; close all  %CANC (?)

global K; global CL_Nominal; global save_figure_flag;

% Carico modello, funzioni peso, ecc
olp_hdd;

% %WC
% %OS %Overshoot
% %AMP
% ST = ST/1000; % Settling Time (NB: converte s in ms)
% %GM  % Gain Margin
% %PM  % Phase Margin
% 
% G_mia=G(2);
% %G_mia=tf(G(2));
% 
% G_mia.InputName = 'u';
% G_mia.OutputName = 'y';
% % G.InputName = {'L','V'};
% % G.OutputName = {'yD','yB'};
% 
% % DM = tunableGain('Decoupler',diag([1 -1]));
% % DM.Gain.Free = [false true;true false];
% 
% %Sum = sumblk('e = r-y');
% 
% % Use the tunablePID object to parameterize the two PI controllers:
% C0 = tunablePID('C0','pid');
% %C0 = tunablePID('C0','pid');
% pid(C0)
% get(C0)
% % PI_L = tunablePID('PI_L','pi');
% % PI_V = tunablePID('PI_V','pi');
% 
% X1 = AnalysisPoint('X1');
% X2 = AnalysisPoint('X2');
% %X3 = AnalysisPoint('X3');
% 
% %CL = feedback(X2*G*C0,X1);
% %CL.InputName = 'r';
% %CL.OutputName = 'y';
% 
% InnerLoop = X2*G_mia*C0;    %in realt? non ? un anello, bens? una connesione serie
% CL0 = feedback(InnerLoop,X1);
% CL0.InputName = 'r';
% CL0.OutputName = 'y';
% 
% % Next construct a model C0 of the controller  in Figure 2.
% %C0 = PI0;
% %%C0 = PI0 * [1, -1];
% %C0 = blkdiag(PI_L,PI_V);
% %C0 = blkdiag(PI_L,PI_V) * DM * [eye(2) -eye(2)];
% 
% %%ss(C0)
% %%get(C0)
% 
% 
% 
% % Note: I/O names should be consistent with those of G
% %C0.InputName = 'e = r - y'; % Dsp/Bsp sono i ref, yD,yB sono le uscite in retroazione
% %%C0.InputName = {'r','y'}; % Dsp/Bsp sono i ref, yD,yB sono le uscite in retroazione
% %%C0.OutputName = 'u'; % controlli
% % C0.InputName = {'Dsp','Bsp','yD','yB'}; % Dsp/Bsp sono i ref, yD,yB sono le uscite in retroazione
% % C0.OutputName = {'L','V'}; % controlli
% 
% % Now tune the controller parameters with looptune as done previously.
% % Crossover frequency
% wc = WC; % 0.5; 
% 
% % Overshoot and disturbance rejection requirements
% %%OS = TuningGoal.Overshoot('r','y',OS); %15
% %%DR = TuningGoal.StepRejection('r','y',AMP,ST); %4, 20 min
% % OS = TuningGoal.Overshoot({'Dsp','Bsp'},{'yD','yB'},OS); %15
% % DR = TuningGoal.StepRejection({'L','V'},{'yD','yB'},AMP,ST); %4, 20 min
% 
% 
% % SOFT requirements
% %Specify tuning requirements for reference tracking and disturbance rejection.
% Rtrack = TuningGoal.Tracking('r','y',ST,0.05);   %the signal at 'y' tracks the signal at 'r' with a response time of ST ms and a tracking error of 1%.
% %Rtrack = TuningGoal.Tracking('r','y',0.001,0.01);   %the signal at 'y' tracks the signal at 'r' with a response time of 5 seconds and a tracking error of 1%.
% Rreject = TuningGoal.Gain('X2','y',1.9);    %limits the gain from the implicit input associated with the AnalysisPoint block X2 to the output, 'y'. (See AnalysisPoint.) Limiting this gain to a value less than 1 ensures that a disturbance injected at X2 is suppressed at the output.
% %Rreject = TuningGoal.Gain('X2','y',0.9);    %limits the gain from the implicit input associated with the AnalysisPoint block X2 to the output, 'y'. (See AnalysisPoint.) Limiting this gain to a value less than 1 ensures that a disturbance injected at X2 is suppressed at the output.
% 
% % HARD requirements
% %Specify tuning requirements for the gain and phase margins.
% RmargOut = TuningGoal.Margins('X1',5,40);
% RmargIn = TuningGoal.Margins('X2',5,40);
% RmargIn.Openings = 'X1';
% % RmargOut = TuningGoal.Margins('X1',18,60);
% % RmargIn = TuningGoal.Margins('X2',18,60);
% % RmargIn.Openings = 'X1';
% 
% ReqOvershoot = TuningGoal.Overshoot('r','y',OS);
% %ReqGain = TuningGoal.Gain('X3','y',1.2)
% 
% % Tune the control system to meet the soft requirements of tracking and disturbance rejection, subject to the hard constraints of the stability margins.
% %SoftReqs = [];
% %SoftReqs = [Rtrack];
% %SoftReqs = [Rreject];
% SoftReqs = [];
% %SoftReqs = [RmargIn,RmargOut];
% %HardReqs = [ReqOvershoot];
% %HardReqs = [RmargIn,RmargOut, ReqOvershoot, ReqGain];
% %HardReqs = [RmargIn,RmargOut, ReqOvershoot];
% HardReqs = [Rtrack,RmargIn,RmargOut];
% %HardReqs = [Rtrack];
% %[CL_final,fSoft,gHard] = systune(prescale(CL0),SoftReqs,HardReqs); %CANC
% [CL_final,fSoft,gHard] = systune(CL0,SoftReqs,HardReqs);
% 
% 
% % Validate the tuned control system against the stability margin requirements.
% %figure;
% %viewGoal(HardReqs,CL_final)   % richiede Matlab 2018b
% 
% % Examine whether the tuned control system meets the tracking requirement by examining the step response from 'r' to 'y'.
% figure;
% stepplot(CL_final,1)
% 
% % Examine the tracking requirement
% %figure;
% %viewGoal(Rtrack,CL_final)     % richiede Matlab 2018b
% 
% 
% 
% 
% 
% 
% % Tune controller gains
% %[~,C,~,Info] = looptune(G,C0,wc);
% %[~,C] = looptune(G,C0,wc,OS);
% %%[~,C] = looptune(G,C0,wc,OS,DR);
% 
% 
% C = CL_final;
% ss(C)
% get(C)
% 
% % Final: Peak gain = 0.997, Iterations = 81
% % Achieved target gain value TargetGain=1.
% % To validate the design, close the loop with the tuned compensator C and simulate the step responses for setpoint tracking and disturbance rejection.
% 
% 
% Tcl = CL_final;
% %Tcl = CL0;
% %Tcl = connect(CL0,CL_final,'r','y');
% %Tcl = connect(G,C0,'r','y');
% %Tcl = connect(G,C,Sum,'r','y');
% %Tcl = connect(G,C,{'Dsp','Bsp','L','V'},{'yD','yB'});
% 
% figure('Position',[0,0,700,350])
% 
% % subplot(121)
% Ttrack = Tcl;
% %Ttrack = Tcl(:,[1 2]);
% step(Ttrack,1), grid, title('Setpoint tracking')
% %step(Ttrack,70), grid, title('Setpoint tracking')
% 
% % subplot(122)
% % Treject = Tcl(:,[3 4]);
% % Treject.InputName = {'dL','dV'};
% % step(Treject, 180), grid, title('Disturbance rejection')
% 
% 
% % K=ss(C);
% % %K=ss(C);
% % %K=ss(PI0);
% % ak = K.A; bk = K.B; ck = K.C; dk = K.D;
% % ak
% % bk
% % ck
% % dk
% % %bk(:,2) = -bk(:,2); %inutile
% % %dk(:,2) = -dk(:,2); %inutile
% % %bk_new=bk(:,1)  %forse non va
% % %dk_new=dk(:,1)  %forse non va
% % %bk(:,3:4) = -bk(:,3:4);    %?
% % %dk(:,3:4) = -dk(:,3:4);    %?
% % % K = pck(ak,bk,ck,dk);
% % % 
% % % showTunable(K)
% % %K=tf(ss(ak,bk_new,ck,dk_new));
% % K=tf(ss(ak,bk,ck,dk));
% K = ss(CL_final.Blocks.C0); %CANC??
% 
% %clp_poles = pole(CL_mu)
% clp_ic = lft(sys_ic,K,1,1);
% clp_poles = pole(clp_ic)
% 
% %tf(CL_final.Blocks.C0)
% %zpk(CL_final.Blocks.C0)
% 
% figure,
% %sigma(K, logspace(-3,3,100)), grid
% bodemag(K,logspace(1,7,500)), grid
% %xlabel('Frequency (rad/min)'), ylabel('Magnitude')
% title('Frequency response of the controller')
% 
% %ss(K)     %CANC
% %get(K)     % CANC
% %tf(K)
% 
% %figure
% %bode(feedback(G_mia*K,1))
% %[Gm , Pm , wg , wp ] = margin(feedback(G_mia*K,1))
% 
% %margin(feedback(G_mia*K,1))
% 
% %clp_poles = pole(CL_final)
% %clp_poles_bo = pole(feedback(G_mia*K,1))
% 
% 
% %isstable(feedback(G_mia*K,1))
% %isstable(CL_final)


% I valori sono stati ricavati eseguendo manualmente il comando Martlab pidTuner
Kp = 0.009603758578229;
Ki = 0.543465518631468;
Kd = 9.370822918985176e-06;
Tf = 1.387549537218259e-05;
Ts = 0;

K_PID = ss(pid(Kp,Ki,Kd,Tf,Ts));
K = K_PID;


figure,
%sigma(K_lsh, logspace(-3,3,100)), grid
bodemag(K_PID,logspace(1,7,500)), grid
%xlabel('Frequency (rad/min)'), ylabel('Magnitude')
title('Frequency response of the controller')
if (save_figure_flag == 1)
    fig_name = ['K_PID__Frequency_response_of_the_controller.png'];
    saveas(gcf, fig_name)
end


clp_ic = lft(sys_ic,K,1,1);
clp_poles = pole(clp_ic)

%end









%%
% function [K,Tcl] = TraditionalDecouplPID( G, WC, OS, AMP, ST )
% 
% G=tf(G);
% 
% G.InputName = {'L','V'};
% G.OutputName = {'yD','yB'};
% 
% DM = tunableGain('Decoupler',diag([1 -1]));
% DM.Gain.Free = [false true;true false];
% 
% 
% % Use the tunablePID object to parameterize the two PI controllers:
% PI_L = tunablePID('PI_L','pi');
% PI_V = tunablePID('PI_V','pi');
% 
% % Next construct a model C0 of the controller  in Figure 2.
% C0 = blkdiag(PI_L,PI_V) * DM * [eye(2) -eye(2)];
% 
% % Note: I/O names should be consistent with those of G
% C0.InputName = {'Dsp','Bsp','yD','yB'}; % Dsp/Bsp sono i ref, yD,yB sono le uscite in retroazione
% C0.OutputName = {'L','V'}; % controlli
% 
% % Now tune the controller parameters with looptune as done previously.
% % Crossover frequency
% wc = WC; % 0.5; 
% 
% % Overshoot and disturbance rejection requirements
% OS = TuningGoal.Overshoot({'Dsp','Bsp'},{'yD','yB'},OS); %15
% DR = TuningGoal.StepRejection({'L','V'},{'yD','yB'},AMP,ST); %4, 20 min
% 
% % Tune controller gains
% [~,C] = looptune(G,C0,wc,OS,DR);
% 
% % Final: Peak gain = 0.997, Iterations = 81
% % Achieved target gain value TargetGain=1.
% % To validate the design, close the loop with the tuned compensator C and simulate the step responses for setpoint tracking and disturbance rejection.
% 
% Tcl = connect(G,C,{'Dsp','Bsp','L','V'},{'yD','yB'});
% 
% figure('Position',[0,0,700,350])
% 
% % subplot(121)
% Ttrack = Tcl(:,[1 2]);
% step(Ttrack,70), grid, title('Setpoint tracking')
% 
% % subplot(122)
% % Treject = Tcl(:,[3 4]);
% % Treject.InputName = {'dL','dV'};
% % step(Treject, 180), grid, title('Disturbance rejection')
% 
% K=ss(C);
% ak = K.A; bk = K.B; ck = K.C; dk = K.D;
% bk(:,3:4) = -bk(:,3:4);
% dk(:,3:4) = -dk(:,3:4);
% % K = pck(ak,bk,ck,dk);
% % 
% % showTunable(K)
% K=tf(ss(ak,bk,ck,dk));
% 
% figure,
% sigma(K, logspace(-3,3,100)), grid
% xlabel('Frequency (rad/min)'), ylabel('Magnitude')
% title('Frequency responses of the controller')
% 
% end