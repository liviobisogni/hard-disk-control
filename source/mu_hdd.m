function [stabmarg, perfmarg] = mu_hdd(sys_ic, K, controller_name)

clc; close all

global save_figure_flag;

% Robustness analysis of the Hard Disk Drive Servo System
%
%clp_ic = lft(sys_ic,K,1,1);
%w = logspace(3,5,100);
%%w = logspace(2,5,100);
%%clp_g = ufrd(clp_ic,w);

clp_ic = lft(sys_ic,K,1,1);
w = logspace(2,5,100);
clp_g = ufrd(clp_ic,w);


% robust stability analysis

% [M,Delta,blkstruct] = lftdata(clp_ic);
% Mfrd = frd(M,w);
% rbnds = mussv(Mfrd(1:15,1:15),blkstruct);

opt = robopt('Display','on');
% opt = robopt('Display','on','Mussv',rbnds);
[stabmarg,destabu,report,info] = robuststab(clp_g,opt);
stabmarg
report
figure(1)
semilogx(info.MussvBnds(1,1),'r-',info.MussvBnds(1,2),'b--')
grid
title('Robust stability')
xlabel('Frequency (rad/s)')
ylabel('mu')
legend('\mu-upper bound','\mu-lower bound','Location','SouthWest') %legend('\mu-upper bound','\mu-lower bound',3)
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Robust_Stability.png'];
    saveas(gcf, fig_name)
end

% nominal performance
w = logspace(2,5,100);
figure(2)
sv = sigma(clp_ic.Nominal,w);
sys_frd = frd(sv(1,:),w);
semilogx(sys_frd,'r-')
grid
title('Nominal performance')
xlabel('Frequency (rad/s)')
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Nominal_Performance.png'];
    saveas(gcf, fig_name)
end



% robust performance
%w = logspace(3,5,100);

%w = logspace(2,5,100);
rp = ucomplexm('rp',zeros(3,2));
%
% mu-controller
%clp_ic = lft(sys_ic,K,1,1);
clp2 = lft(clp_ic,rp);
[M,Delta,blkstruct] = lftdata(clp2);
Mfrd = frd(M,w);
bnd = mussv(Mfrd(1:17,1:18),blkstruct);



%w = logspace(2,5,100);
opt = robopt('Display','on');
% opt = robopt('Display','on','Mussv',bnd);
[perfmarg,perfmargunc,report,info] = robustperf(clp_g,opt);
perfmarg
report
figure(3)
semilogx(info.MussvBnds(1,1),'r-',info.MussvBnds(1,2),'b--')
grid
title('Robust performance')
xlabel('Frequency (rad/s)')
ylabel('mu')
legend('\mu-upper bound','\mu-lower bound','Location','SouthWest')   %legend('\mu-upper bound','\mu-lower bound',3)
if (save_figure_flag == 1)
    fig_name = [controller_name, '__Robust_Performance.png'];
    saveas(gcf, fig_name)
end
end