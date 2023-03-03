function [] = verifiche(controller_name)

% Verifiche

%          * Peak closed-loop gain < 4 dB
%          * Open-loop gain > 20 dB at 100 Hz 
%          * Steady-state error < 0.1 ?m (micrometro)
%          * Settling time < 1.5 ms
%          * Closed-loop bandwidth > 1000 Hz
%          * Gain margin > 5 dB
%          * Phase margin > 40 deg
%          * u(t) < 1.2 V

close all; clc
global K; global CL_Nominal; global save_figure_flag;

sim_hdd
clp = lft(sim_ic,K,1,1);    % Input Name: ref, dist, noise; Output Name: G(1), control

% closed-loop Bode plot                                 % ref G(1)
ref_loop = clp(1,1);

% output sensitivity to disturbance                     % dist G(1)
sen_loop = clp(1,2);

% output sensitivity to noise                           % noise G(1)
noise_loop = clp(1,3);

% sensitivity of the control effort to the reference    % ref u
ref_u = clp(2,1);

% sensitivity of the control effort to the disturbance  % dist u
dist_u = clp(2,2);

% sensitivity of the control effort to the noise        % noise u
noise_u = clp(2,3);



%          * Peak closed-loop gain < 4 dB
%ref_loop = clp(1,1);   %CANC
w = logspace(2,5,400);
figure(1)
bodemag(ref_loop.Nominal,'r-',ref_loop,'b--',w), grid % closed-loop Bode plot
%bode(ref_loop.Nominal,'r-',ref_loop,'b--',w), grid % closed-loop Bode plot
title('Closed-loop Bode plot')
legend('Nominal system','Random samples','Location','SouthWest') %legend('Nominal system','Random samples',3)
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Closed-loop_Bode_plot.png'];
    saveas(gcf, fig_name)
end
[gpeak,fpeak] = getPeakGain(ref_loop.Nominal);
Peak_closed_loop_gain = mag2db(gpeak);   %dB
if (Peak_closed_loop_gain < 4)
    temp = ['Peak closed-loop gain = ', num2str(Peak_closed_loop_gain), ' dB < 4 dB: specifica soddisfatta'];
else
    temp = ['Peak closed-loop gain = ', num2str(Peak_closed_loop_gain), ' dB >= 4 dB: specifica NON soddisfatta'];
end
disp(temp)
%abs(freqresp(ref_loop,4.7796e+03))



%          * Open-loop gain > 20 dB at 100 Hz 
L = G(1,2)*K;
w = logspace(0,5,400);
figure(2)
bode(L.Nominal,'r-',L,'b--',w), grid % open-loop frequency response
bodemag(L.Nominal,'r-',L,'b--',w), grid % open-loop frequency response
title('Open-loop Bode plot')
legend('Nominal system','Random samples','Location','SouthWest')    %legend('Nominal system','Random samples',3)
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Open-loop_Bode_plot.png'];
    saveas(gcf, fig_name)
end
%freqresp(L,4.7796e+03)
Open_loop_gain_at_100Hz = mag2db(abs(freqresp(L,100/(2*pi))));  %dB
if (Open_loop_gain_at_100Hz > 20)
    temp = ['Open-loop gain(100 Hz) = ', num2str(Open_loop_gain_at_100Hz), ' dB > 20 dB: specifica soddisfatta'];
else
    temp = ['Open-loop gain(100 Hz) = ', num2str(Open_loop_gain_at_100Hz), ' dB <= 20 dB: specifica NON soddisfatta'];
end
disp(temp)



stepinfo_temp = stepinfo(ref_loop); %
%          * Steady-state error < 0.1 micro m
%stepinfo_temp = stepinfo(clp(1,1)) %CANC
SP=1; %input value, if you put 1 then is the same as step(sys)
figure(3)
%step(SP*ref_loop.NominalValue,'k--')
%step(SP*ref_loop.NominalValue,'k--', 'LineWidth', 2)
%hold on

step(SP*ref_loop.NominalValue, 'k--')
hold on
step(SP*ref_loop, 'c-')
step(SP*ref_loop.NominalValue, 'k--')


% [y,t] = step(SP*ref_loop);
% [y_nom,t_nom] = step(SP*ref_loop.NominalValue);
% plot(t, y_nom, 'k--', 'LineWidth', 2)
% hold on
% plot(t, y, 'c-')
% plot(t, y_nom, 'k--', 'LineWidth', 2)

%title('Control action, due to reference')
legend('Nominal system','Random samples','Location','SouthEast') %legend('Nominal system','Random samples',3)
%xlabel('Time (secs)')
%ylabel('u (V)')
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Step_response.png'];
    saveas(gcf, fig_name)
end
[y,t]=step(SP*ref_loop,0.001); %get the response of the system to a step with amplitude SP
sserror_1=abs(SP-y(end));  %!!!      %get the steady state error



T=SP*ref_loop;
format long g;
ss_values = dcgain(T);  % Steady-state values
sserror_2 = abs(SP-ss_values);  %!!!
%disp(ss_values);
if (sserror_2 < 1e-7)
    temp = ['Steady-state error = ', num2str(sserror_2), ' m < 1e-7 m: specifica soddisfatta'];
else
    temp = ['Steady-state error = ', num2str(sserror_2), ' m >= 1e-7 m: specifica NON soddisfatta'];
end
disp(temp)



%          * Settling time < 1.5 ms
stepinfo_temp.SettlingTime;
if (stepinfo_temp.SettlingTime < 1.5e-3)
    temp = ['Settling time = ', num2str(stepinfo_temp.SettlingTime), ' s < 0.0015 s: specifica soddisfatta'];
else
    temp = ['Settling time = ', num2str(stepinfo_temp.SettlingTime), ' s >= 0.0015 s: specifica NON soddisfatta'];
end
disp(temp)



%          * Closed-loop bandwidth > 1000 Hz
%%MANCA!!!
fb = bandwidth(ref_loop);    % rad/s
fb_Hz = fb*2*pi;            % Hz
if (fb_Hz > 1000)
    temp = ['Closed-loop bandwidth = ', num2str(fb_Hz), ' Hz > 1000 Hz: specifica soddisfatta'];
else
    temp = ['Closed-loop bandwidth = ', num2str(fb_Hz), ' Hz <= 1000 Hz: specifica NON soddisfatta'];
end
disp(temp)


figure(4)
margin(G(2).NominalValue*K)
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Gain_and_Phase_margins.png'];
    saveas(gcf, fig_name)
end
[Gm , Pm , wg , wp ] = margin(G(2).NominalValue*K);    %
Gm_dB = mag2db(Gm);
%          * Gain margin > 5 dB
if (Gm_dB > 5)
    temp = ['Gain margin = ', num2str(Gm_dB), ' dB > 5 dB: specifica soddisfatta'];
else
    temp = ['Gain margin = ', num2str(Gm_dB), ' dB <= 5 dB: specifica NON soddisfatta'];
end
disp(temp)

%          * Phase margin > 40 deg
if (Pm > 40)
    temp = ['Phase margin = ', num2str(Pm), ' deg > 40 deg: specifica soddisfatta'];
else
    temp = ['Phase margin = ', num2str(Pm), ' deg <= 40 deg: specifica NON soddisfatta'];
end
disp(temp)



%          * u(t) < 1.2 V
% response to the reference
r = 1.2;                                % 1 track --> 1.2 V 
ti = 0.000001;                          % time increment
tfin = 0.005;                           % final time value
time = 0:ti:tfin;
nstep = size(time,2);
ref(1:nstep) = r;
dist(1:nstep) = 0.0;
noise(1:nstep) = 0.0;
nsample = 30;
clp_30 = usample(clp,nsample);

[y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
yr = y(:,1);
ur_nominal = y(:,2);
%ym = yr/r;
%err_nominal  = ref' - yr;                        % PES in volts
figure(5)
plot(t,ur_nominal,'k--', 'LineWidth', 2)
hold on
for i = 1:nsample
    [y,t] = lsim(clp_30(1:2,1:3,i),[ref',dist',noise'],time);
%     yr = y(:,1);
%     err  = ref' - yr;                   % PES in volts
%     errm = err/r;                       % PES in tracks
%     figure(5)
%     plot(t,errm,'b-')
%     hold on
    ur = y(:,2);
    figure(5)
    plot(t,ur,'r-')
    axis([0 0.0005 -0.6 1.4])
    hold on
end
% figure(5)
% grid
% title('Error transient response to reference')
% xlabel('Time (secs)')
% ylabel('Position Error Signal (tracks)')
% hold off
% figure(5)
% grid
% title('Closed-loop transient response')
% xlabel('Time (secs)')
% ylabel('Position Error Signal (tracks)')
% hold off
[y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
yr = y(:,1);
ur_nominal = y(:,2);
%ym = yr/r;
%err_nominal  = ref' - yr;                        % PES in volts
figure(5)
plot(t,ur_nominal,'k--', 'LineWidth', 2)
grid
title('Control action, due to reference')
legend('Nominal system','Random samples','Location','SouthEast') %legend('Nominal system','Random samples',3)
xlabel('Time (secs)')
ylabel('u (V)')
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Control_action_due_to_reference.png'];
    saveas(gcf, fig_name)
end
hold off
clear ref, clear dist, clear noise

max_ur = max(ur);

if (max_ur < 1.2)
    temp = ['max u_r(t) = ', num2str(max_ur), ' V < 1.2 V: specifica soddisfatta'];
else
    temp = ['max u_r(t) = ', num2str(max_ur), ' V >= 1.2 V: specifica NON soddisfatta'];
end
disp(temp)
disp('Fine verifiche')      %CANC





% response to the disturbance
r = 1.2;                                % 1 track --> 1.2 V 
ti = 0.000001;                           % time increment
tfin = 0.005;                           % final time value
time = 0:ti:tfin;
nstep = size(time,2);
ref(1:nstep) = 0.0;
dist(1:nstep) = 0.0005;                 % Td = 0.0005 N.m
noise(1:nstep) = 0.0;
nsample = 30;
clp_30 = usample(clp,nsample);

[y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
yd = y(:,1);
%ud = y(:,2);
%ym = yd/r;
err_nominal  = ref' - yd;                           % PES in volts
errm_nominal = err_nominal/r;                       % PES in tracks
figure(6)
plot(t,errm_nominal,'k--', 'LineWidth', 2)
hold on
for i = 1:nsample
    [y,t] = lsim(clp_30(1:2,1:3,i),[ref',dist',noise'],time);
    yd = y(:,1);
    err  = ref' - yd;                   % PES in volts
    errm = err/r;                       % PES in tracks
    figure(6)
    plot(t,errm,'c-')
     hold on
%     ud = y(:,2);
%     figure(6)
%     plot(t,ud,'r-')
%     axis([0 0.0005 -0.02 0.01])
%     hold on
end

[y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
yd = y(:,1);
%ud = y(:,2);
%ym = yd/r;
err_nominal  = ref' - yd;                           % PES in volts
errm_nominal = err_nominal/r;                       % PES in tracks
figure(6)
plot(t,errm_nominal,'k--', 'LineWidth', 2)

grid
title('Error transient response to disturbance')
legend('Nominal system','Random samples','Location','SouthEast') %legend('Nominal system','Random samples',3)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
% text(xlim(2)-1e-3,ylim(2)-1e-2,'ref = 0')
% text(xlim(2)-1e-3,ylim(2)-1.7e-2,'dist = 0.0005 N*m')
% text(xlim(2)-1e-3,ylim(2)-2.4e-2,'noise = 0')
% text(xlim(2)-1e-3,ylim(2)-ylim(2)/8,'ref = 0')
% text(xlim(2)-1e-3,ylim(2)-ylim(2)/2.5,'dist = 0.0005 N*m')
% text(xlim(2)-1e-3,ylim(2)/3,'noise = 0')
text(xlim(2)-1e-3,ylim(2)-ylim(2)/8,'ref = 0')
text(xlim(2)-1e-3,2*ylim(2)/3-ylim(2)/8,'dist = 0.0005 N*m')
text(xlim(2)-1e-3,ylim(2)/3-ylim(2)/8,'noise = 0')
xlabel('Time (secs)')
ylabel('Position Error Signal (tracks)')
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Error_transient_response_to_disturbance.png'];
    saveas(gcf, fig_name)
end
hold off
% figure(6)
% grid
% title('Transient response to disturbance')
% xlabel('Time (secs)')
% ylabel('Position Error Signal (tracks)')
% hold off
% figure(6)
% grid
% title('Control action due to disturbance')
% xlabel('Time (secs)')
% ylabel('u (V)')
% hold off
% clear ref, clear dist, clear noise





% response to noise
rand('state',0)
r = 1.2;                                % 1 track --> 1.2 V 
ti = 0.000001;                          % time increment
tfin = 0.005;                           % final time value
time = 0:ti:tfin;
nstep = size(time,2);
ref(1:nstep) = 0.0;
dist(1:nstep) = 0.0;                 
noise(1:nstep) = 2.0*(rand(1,nstep)-0.5*ones(1,nstep));

nsample = 30;
clp_30 = usample(clp,nsample);

[y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
yd = y(:,1);
%ud = y(:,2);
%ym = yd/r;
err_nominal  = ref' - yd;                           % PES in volts
errm_nominal = err_nominal/r;                       % PES in tracks
figure(7)
plot(t,errm_nominal,'k--', 'LineWidth', 2)
hold on
for i = 1:nsample
    [y,t] = lsim(clp_30(1:2,1:3,i),[ref',dist',noise'],time);
    yr = y(:,1);
    err  = ref' - yr;                   % PES in volts
    errm = err/r;                       % PES in tracks
    figure(7)
    plot(t,errm,'c-')
     hold on
%     ur = y(:,2);
%     figure(7)
%     plot(t,ur,'r-')
%     axis([0 0.0005 -0.6 1.4])
%     hold on
end

[y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
yd = y(:,1);
%ud = y(:,2);
%ym = yd/r;
err_nominal  = ref' - yd;                           % PES in volts
errm_nominal = err_nominal/r;                       % PES in tracks
figure(7)
plot(t,errm_nominal,'k--', 'LineWidth', 2)
hold on

% [y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
% yd = y(:,1);
% ud = y(:,2);
% ym = yd/r;
% err_nominal  = ref' - yd;                       % PES in volts
% %errm = err/r;                           % PES in tracks


grid
title('Error transient response to noise')
legend('Nominal system','Random samples','Location','NorthWest') %legend('Nominal system','Random samples',3)
%text('ref = 0', 'dist = 0');
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
% text(xlim(2)-1e-3,ylim(2)-1e-3,'ref = 0')
% text(xlim(2)-1e-3,ylim(2)-1.7e-3,'dist = 0')
% text(xlim(2)-1e-3,ylim(2)-2.4e-3,'-1 < noise < 1')
text(xlim(2)-1e-3,ylim(2)-ylim(2)/8,'ref = 0')
text(xlim(2)-1e-3,ylim(2)-ylim(2)/4.7,'dist = 0')
text(xlim(2)-1e-3,ylim(2)-ylim(2)/3.3,'-1 < noise < 1')

% text(xlim(2)-1e-3,ylim(2)-ylim(2)/8,'ref = 0')
% text(xlim(2)-1e-3,2*ylim(2)/3-ylim(2)/8,'dist = 0.0005 N*m')
% text(xlim(2)-1e-3,ylim(2)/3-ylim(2)/8,'noise = 0')
xlabel('Time (secs)')
ylabel('Position Error Signal (tracks)')
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Error_transient_response_to_noise.png'];
    saveas(gcf, fig_name)
end

end
% figure(7)
% grid
% title('Transient response to disturbance')
% xlabel('Time (secs)')
% ylabel('Position Error Signal (tracks)')
% hold off
% figure(7)
% grid
% title('Control action due to disturbance')
% xlabel('Time (secs)')
% ylabel('u (V)')
% hold off
% clear ref, clear dist, clear noise




% [y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
% yn = y(:,1);
% un = y(:,2);
% ym = yn/r;
% err  = ref' - yn;                       % PES in volts
% errm = err/r;                           % PES in tracks
% figure(5)
% plot(t,errm,'c-')
% grid
% title('Transient response to noise')
% xlabel('Time (secs)')
% ylabel('Position Error Signal (tracks)')
% disp(['noise error: ' num2str(100*(norm(ym,inf))) '%'])
% figure(6)
% plot(t,un,'r-')
% grid
% title('Control action due to noise')
% xlabel('Time (secs)')
% ylabel('u (V)')
% clear ref, clear dist, clear noise








% % output sensitivity to disturbance                     % dist G(1)
% figure(6)
% step(sen_loop)
% grid
% title('output sensitivity to disturbance')
% xlabel('Time (secs)')
% ylabel(' (tracks)')
% hold off
% 
% 
% % output sensitivity to noise                           % noise G(1)
% figure(7)
% step(noise_loop)
% grid
% title('output sensitivity to noise')
% xlabel('Time (secs)')
% ylabel(' (tracks)')
% hold off





% 
% %
% % response to the reference
% r = 1.2;                                % 1 track --> 1.2 V 
% ti = 0.000001;                          % time increment
% tfin = 0.005;                           % final time value
% time = 0:ti:tfin;
% nstep = size(time,2);
% ref(1:nstep) = r;
% dist(1:nstep) = 0.0;
% noise(1:nstep) = 0.0;
% [y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
% yr = y(:,1);
% ur = y(:,2);
% ym = yr/r;
% err  = ref' - yr;                        % PES in volts
% errm = err/r;                            % PES in tracks
% figure(9)
% plot(t,errm,'c-')
% grid
% title('Closed-loop transient response')
% xlabel('Time (secs)')
% ylabel('Position Error Signal (tracks)')
% disp(['overshoot: ' num2str(100*(norm(ym,inf)-1)) '%'])
% figure(10)
% plot(t,y(:,2),'r-')
% axis([0 0.0005 -0.6 1.4])
% grid
% title('Control action, due to reference')
% xlabel('Time (secs)')
% ylabel('u (V)')
% disp(['u_max: ' num2str(norm(ur,inf)) 'V'])
% clear ref, clear dist, clear noise
% %
% % response to the disturbance
% r = 1.2;                                % 1 track --> 1.2 V 
% ti = 0.000001;                          % time increment
% tfin = 0.005;                           % final time value
% time = 0:ti:tfin;
% nstep = size(time,2);
% ref(1:nstep) = 0.0;
% dist(1:nstep) = 0.0005;                 % Td = 0.0005 N.m
% noise(1:nstep) = 0.0;
% [y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
% yd = y(:,1);
% ud = y(:,2);
% ym = yd/r;
% err  = ref' - yd;                       % PES in volts
% errm = err/r;                           % PES in tracks
% figure(11)
% plot(t,errm,'c-')
% grid
% title('Transient response to disturbance')
% xlabel('Time (secs)')
% ylabel('Position Error Signal (tracks)')
% disp(['dist. error: ' num2str(100*(norm(ym,inf))) '%'])
% figure(12)
% plot(t,ud,'r-')
% grid
% title('Control action due to disturbance')
% xlabel('Time (secs)')
% ylabel('u (V)')
% clear ref, clear dist, clear noise
% %
% % response to noise
% rand('state',0)
% r = 1.2;                                % 1 track --> 1.2 V 
% ti = 0.000001;                          % time increment
% tfin = 0.005;                           % final time value
% time = 0:ti:tfin;
% nstep = size(time,2);
% ref(1:nstep) = 0.0;
% dist(1:nstep) = 0.0;                 
% noise(1:nstep) = 2.0*(rand(1,nstep)-0.5*ones(1,nstep));
% [y,t] = lsim(clp.Nominal(1:2,1:3),[ref',dist',noise'],time);
% yn = y(:,1);
% un = y(:,2);
% ym = yn/r;
% err  = ref' - yn;                       % PES in volts
% errm = err/r;                           % PES in tracks
% figure(13)
% plot(t,errm,'c-')
% grid
% title('Transient response to noise')
% xlabel('Time (secs)')
% ylabel('Position Error Signal (tracks)')
% disp(['noise error: ' num2str(100*(norm(ym,inf))) '%'])
% figure(14)
% plot(t,un,'r-')
% grid
% title('Control action due to noise')
% xlabel('Time (secs)')
% ylabel('u (V)')
% clear ref, clear dist, clear noise