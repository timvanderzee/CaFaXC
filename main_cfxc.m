clear all; close all; clc
saveFig = 1;

%% Code block: set Path.
%%% set the codefolder variable which is all the pathing that needs to be done.
% make sure that you're in the same folder as the main_cfxc.m file
[codefolder, name, ext] = fileparts(which('main_cfxc.m')); % folder contains the present m-file
% which we expect to be the current working directory.

savedir = fullfile(codefolder,'manuscript figures','model output');

% add sub-folders to path that contain data:mai
% 'kinematics_kinetics.mat' and 'muscle_excitations.mat'
addpath(fullfile(codefolder));
addpath(fullfile(codefolder,'model parameters'));
addpath(fullfile(codefolder,'cyclic force study','inputs'));

%% Update parameters
setts.mouse = 0; % false for human

if setts.mouse
    load('mouse_parms.mat')
    parms.exp.phi = 120; % joint angle (full extension = 0)
else
    load('quad_parms.mat','parms','fl','fv');
    parms.exp.phi = 90; % joint angle (full extension = 0)
end

% update
parms.ce.amin = 1e-3; % minimal excitation

parms = cfxc.calc_x0(parms); 

% needs to be stored here, because will be changed in next section
tau_original = parms.ce.tau;

%% Contributions to force development (tetanus)
parms.exp.stim_type = 'constant';
parms.exp.a = 1;
parms.exp.A = 1;

parms.set.optimum = 1;

% simulate isometric contraction
parms.set.no_tendon = 0;
[t1,x1] = ode113(@cfxc.sim_muscle, [0 .5], [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end)], parms.set.odeopt, parms);

% to isolate contraction from crossbridge, remove the tendon
parms.set.no_tendon = 1;
[t2,x2] = ode113(@cfxc.sim_muscle, [0 .5], [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end)], parms.set.odeopt, parms);

Ca = x1(:,1);
Fa = x1(:,2);
F = x1(:,4)/parms.CB.Xmax(2);
F2 = x2(:,4)/parms.CB.Xmax(2);

Trise = [min(t1(Ca > .9)) min(t1(Fa > .9))  min(t2(F2 > .9)) min(t1(F > .9))];
cont = round(100*[Trise(1) diff(Trise)] / Trise(end),2);

if ishandle(1), close(1); end; figure(1)
plot(t1,[Ca Fa]); hold on; box off
plot(t2, F2);
plot(t1, F);
yline(0.9,'k--')

legend('Calcium','Facilitation','Force without tendon','Force','location','best'); 
legend boxoff
xlabel('Time (s)'); ylabel('State (a.u.)');

%% Cyclic force production (Fig 8)
setts.mouse = 0; % false for human
load('quad_parms.mat','parms','fl','fv');
load('kinematics_kinetics.mat','freqs','tx','Tx');

parms.set.no_tendon = 0;
parms.set.odeopt = odeset('maxstep',1e-3);
parms.set.optimum = 0;

parms.exp.phi = 20;
parms = cfxc.calc_x0(parms);

% parameters
parms.exp.vmtc = 0;
parms.exp.stim_type = 'interp';

parms.type = 'CaFaXC';

if ishandle(2),close(2); end
figure(2)

% load input
load('muscle_excitations','U_opt_smooth');

colors = parula(length(freqs));

% display the target
for f = 1:length(freqs)
  subplot(234)
  plot(tx(f,:), Tx,'k--'); hold on
end

% pre-allocate
tForce = nan(size(freqs));
tExcit = nan(size(freqs));
Fgain = nan(size(freqs));

for f = 1:length(freqs)

    disp(['Freq = ', num2str(freqs(f)), ' Hz'])
    parms.exp.t = tx(f,:);
    parms.exp.U = U_opt_smooth(f,:);

    % simulate
    [t,x] = ode113(@cfxc.sim_muscle, [0 max(parms.exp.t)], [0 parms.exp.x0], parms.set.odeopt, parms);
      
    F = x(:,4)/parms.CB.delta;
    lce = x(:,6);

    W = cumtrapz(-lce,F);
    
    subplot(231);
    plot(parms.exp.t, parms.exp.U,'color',colors(f,:)); hold on;
    
    subplot(232); 
    plot(t,x(:,1),'color',colors(f,:)); hold on;
    
    subplot(233); 
    plot(t,x(:,2),'color',colors(f,:)); hold on; 
    
    subplot(234)
    plot(t,F * parms.mtc.r,'color',colors(f,:)); hold on

    subplot(235);
    plot(t,lce,'color',colors(f,:)); hold on
    
    subplot(236);
    plot(t,W,'color',colors(f,:)); hold on
    
    % low-pass filter excitation
    tup = linspace(0,max(parms.exp.t),1000);
    Uup = interp1(parms.exp.t, parms.exp.U, tup);
    
    Wn = 10 / (.5/mean(diff(tup)));
    [b,a] = butter(2, Wn);
    Ufilt = filtfilt(b,a,Uup); 
    
    subplot(231)
    plot(tup, Ufilt, ':','color', colors(f,:))
        
    tForce(f) = t(F==max(F));
    tExcit(f) = tup(Ufilt==(max(Ufilt)));

    Fgain(f) = 1/max(Uup);
end

Fphase = -(tForce-tExcit).* freqs * 360;

% make nice
titles = {'Excitation', 'Calcium activation', 'Force activation', 'Force','Length','Work'};
ylabels = {'u (0-1)', 'a_{Ca} (0-1)', 'a_F (0-1)', 'Torque (N-m)', 'Length (m)', 'Work (J)'};

for i = 1:6
    subplot(2,3,i);
    box off
    title(titles{i})
    xlabel('Time (s)')
    ylabel(ylabels{i})
end

set(gcf,'units','normalized','position', [0 .3 .6 .4])

%% summary figure
if ishandle(3), close(3); end; figure(3)
subplot(211); semilogx(freqs, Fgain/Fgain(1),'linewidth',2); box off; hold on
ylim([0 1.5])
xlabel('Cyclic contraction frequency (Hz)'); ylabel('Gain (a.u.)'); title('Gain')
xlim([.5 3.1])

subplot(212); semilogx(freqs, Fphase,'linewidth',2); box off; hold on
ylim([-150 0])
xlabel('Cyclic contraction frequency (Hz)'); ylabel('Phase (deg)'); title('Phase')
xlim([.5 3.1])

set(gcf,'units','normalized','position', [.6 .3 .2 .4])

if saveFig
    save(fullfile(savedir,'Fig8','Fig8_lowpass-filter.mat'),'freqs','Fgain','Fphase','parms')
end

%% Force-frequency (Fig 8)
parms.exp.phi = 90;
parms = cfxc.calc_x0(parms);

parms.exp.stim_type = 'u_func';
parms.exp.A = 1;
parms.type = 'CaFaXC';

freqs = logspace(0,3,50);
tmax = 0.3;
parms.exp.tstop = tmax;

colors = parula(length(freqs));

if ishandle(4), close(4); end; figure(4)

Fpeak_freqs = nan(size(freqs));

parms.exp.u_func = @(t, parms) parms.exp.A .* (.5 + .5 * square(2*pi*parms.exp.freq*t, .5 * parms.exp.freq)) .* (t < parms.exp.tstop);
    
for f = 1:length(freqs)
    
    disp(['Freq = ', num2str(freqs(f)), ' Hz'])
        
    parms.exp.freq = freqs(f);
    
     % simulate
    [t,x] = ode113(@cfxc.sim_muscle, [0 tmax], [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end)], parms.set.odeopt, parms);

    Ca = x(:,1);
    Fac = x(:,2);
    lce = x(:,6);
    Fce = x(:,4) / parms.CB.delta;
   
    Fpe = parms.func.fpe(lce, parms);
    lse = parms.exp.lmtc - lce;
    Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;

    W = cumtrapz(-lce,Fce);
    
    figure(4)
    
    subplot(231);
    plot(t, parms.exp.u_func(t,parms),'color',colors(f,:)); hold on;
    
    subplot(232); 
    plot(t,Ca,'color',colors(f,:)); hold on;
    
    subplot(233); 
    plot(t,Fac,'color',colors(f,:)); hold on; 
    
    subplot(234)
    plot(t,Fce,'color',colors(f,:)); hold on
    plot(t,Fse,':','color',colors(f,:));
    plot(t,Fce+Fpe,'--','color',colors(f,:));

    subplot(235);
    plot(t,lce-lce(1),'color',colors(f,:)); hold on
    
    subplot(236);
    plot(t,W,'color',colors(f,:)); hold on
        
    Fpeak_freqs(f) = max(Fse) - min(Fse);
end

% make nice
titles = {'Excitation', 'Calcium activation', 'Force activation', 'Force','Length','Work'};
ylabels = {'u (0-1)', 'a_{Ca} (0-1)', 'a_F (0-1)', 'Force (N)', 'Length (m)', 'Work (J)'};

for i = 1:6
    subplot(2,3,i);
    box off
    title(titles{i})
    xlabel('Time (s)')
    ylabel(ylabels{i})
end

set(gcf,'units','normalized','position', [0 .3 .6 .4])

% summary figure
if ishandle(5), close(5); end; figure(5)
plot(freqs, Fpeak_freqs,'linewidth',2)
axis([0 100 0 max(Fpeak_freqs)])
box off
xlabel('Stimulation rate (Hz)'); ylabel('Peak force (N)')
title('Force - frequency')

set(gcf,'units','normalized','position', [.6 .3 .2 .4])

if saveFig
    save(fullfile(savedir,'Fig8','Fig8_force-frequency.mat'),'freqs','Fpeak_freqs','parms')
end

%% Force-velocity (Fig 8)
load('quad_parms.mat')
% inspired by the experiment from Westing et al. (1990)
% Acta Physiol Scand, 140, 17-22
if ishandle(6), close(6); end; figure(6)
% time the muscle is held isometric before shortening is allowed. not sure
% whether Westing is doing this, other should be set to 0. but other
% experiments may have a hold period before releasing
tiso = 1; % [s]

vall = (-parms.ce.vmaxrel:.5:(parms.ce.vmaxrel/2)) * parms.ce.lceopt;
vall = sort(unique([vall 0])); % make sure to include 0

Fss = nan(size(vall));
colors = parula(length(vall));

parms.set.sim_mtc = 0;
parms.exp.vmtc = 0;
parms = cfxc.calc_x0(parms);

% stim parameters
parms.exp.stim_type = 'max';
parms.set.no_tendon = 0;
parms.set.optimum = 0;
parms.CB.analytical = 1;

parms.set.odeopt = odeset('maxstep',1e-3);
        
k = 0;

% isometric angle (yielding optimum length)
phi_iso = 80; % [deg]

% angle excursion
dphi = 45; % [deg]

for c = 1:3

    if c == 1 % concentric 
        parms.exp.phi = 120;
        vs = vall(vall<0);
    elseif c == 2 % isometric
        parms.exp.phi = phi_iso;
        vs = 0; 
    elseif c == 3 % eccentric
        parms.exp.phi = 40;
        vs = vall(vall>0);
    end
    
    % for updating parameters
    parms.set.sim_mtc = 0;
    parms.exp.vmtc = 0;
    parms = cfxc.calc_x0(parms);

    % for the force-velocity experiment
    parms.type = 'CaFaXC';
    parms.set.sim_mtc = 1; % simulate MTC
    parms.exp.stim_type = 'max';

for j = 1:length(vs)
    
    disp(['Velocity = ', num2str(vs(j)), ' m/s'])
    
    k = k+1;
    
    % isometric phase
    parms.exp.vmtc = 0;     
    [t0,x0] = ode113(@cfxc.sim_muscle, [0 tiso], [1 parms.exp.x0(1) parms.exp.x0(2:end) parms.exp.lmtc], parms.set.odeopt, parms);
 
    % isokinetic phase
    parms.exp.vmtc = vs(j);    
    phidot = parms.exp.vmtc / parms.mtc.r * 180/pi;  
    tmax = min([abs(dphi / phidot), 1.5]);

    [t1,x1] = ode113(@cfxc.sim_muscle, [0 tmax], x0(end,:), parms.set.odeopt, parms);

%     t = [t0; t1+t0(end)];
%     x = [x0; x1];

    t = t1;
    x = x1;

    Ca = x(:,1);
    Fac = x(:,2);
    lce = x(:,6);
    Fce = x(:,4) / parms.CB.delta;
    lmtc = x(:,7);
    
    Fpe = parms.func.fpe(lce, parms);
    lse = lmtc - lce;
    Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;    
    W = cumtrapz(-lce,Fce);
    
    subplot(231);
    plot(t, Ca,'color',colors(k,:)); hold on;
    
    subplot(232); 
    plot(t,Ca,'color',colors(k,:)); hold on;
    
    subplot(233); 
    plot(t,Fac,'color',colors(k,:)); hold on; 
    
    subplot(234)
    plot(t,Fce,'color',colors(k,:)); hold on
    plot(t,Fse,':','color',colors(k,:));
    plot(t,Fce+Fpe,'--','color',colors(k,:));
    
    subplot(235);
    plot(t,lce,'color',colors(k,:)); hold on
    
    subplot(236);
    plot(t,lmtc,'color',colors(k,:)); hold on

    Fss(k) = Fce(end);
         
end
end

% make nice
titles = {'Excitation', 'Calcium activation', 'Force activation', 'Force','CE length','MTC length'};
ylabels = {'u (0-1)', 'a_{Ca} (0-1)', 'a_F (0-1)', 'Force (N)', 'Length (m)', 'Length (m)'};

for i = 1:6
    subplot(2,3,i);
    box off
    title(titles{i})
    xlabel('Time (s)')
    ylabel(ylabels{i})
    axis tight
end

set(gcf,'units','normalized','position', [0 .3 .6 .4])


%% summary figure
if ishandle(7), close(7); end; figure(7);
color = get(gca,'colororder');

% Hill-type force-velocity
parms.Fasymp = 1.5;
Fse = linspace(0,parms.Fasymp,1000);
Fisom = 1; 
a = 1;
Vce = parms.func.fv(a, Fse, Fisom, parms);

plot(vall, Fss,'-','linewidth',2,'color',color(1,:),'markerfacecolor',color(1,:),'markersize',5); hold on; box off;
plot(Vce*parms.ce.lceopt, Fse * parms.ce.Fmax,'--'); 
plot(fv.vHill*parms.ce.lceopt, fv.FCB(:,[1 3])*parms.ce.Fmax,':')
legend('CaFaXC','Hill','CB (Huxley)','CB (DM)','location','best')
legend boxoff

axis([-parms.ce.vmaxrel*parms.ce.lceopt .5*parms.ce.vmaxrel*parms.ce.lceopt 0 2*parms.ce.Fmax])
xlabel('MTC velocity (m/s)'); ylabel('Force (N)'); 

set(gcf,'units','normalized','position', [.6 .3 .2 .4])
title('Force-velocity')

if saveFig
    save(fullfile(savedir,'Fig8','Fig8_force-velocity.mat'),'Vce','Fse','vall','Fss');
end

%% Force-length (Fig 8)
if ishandle(8), close(8); end; figure(8)

% settings
parms.exp.vmtc = 0;
parms.set.sim_mtc = 0;
parms.ce.amin = 0;

phis = 0:5:120;
colors = parula(length(phis));

% pre-allocate
Fse_max = nan(size(phis));
Fce_max = nan(size(phis));
Lmin = nan(size(phis));
Lmax = nan(size(phis));

for i = 1:length(phis)
    
    disp(['Joint angle = ', num2str(phis(i)), ' deg'])
    
    % update
    parms.exp.phi = phis(i);
    parms = cfxc.calc_x0(parms);
    
    % simulate
    parms.exp.stim_type = 'max';
    parms.type = 'CaFaXC';
    [t,x] = ode113(@cfxc.sim_muscle, [0 .5], [1 parms.exp.x0(1) parms.exp.x0(2:end)], parms.set.odeopt, parms);

    Ca = x(:,1);
    Fac = x(:,2);
    lce = x(:,6);
    Fce = x(:,4) / parms.CB.delta;
    
    Fpe = parms.func.fpe(lce, parms);
    lse = parms.exp.lmtc - lce;
    Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;    
    W = cumtrapz(-lce,Fce);
    
    subplot(231);
    plot(t, Ca,'color',colors(i,:)); hold on;
    
    subplot(232); 
    plot(t,Ca,'color',colors(i,:)); hold on;
    
    subplot(233); 
    plot(t,Fac,'color',colors(i,:)); hold on; 
    
    subplot(234)
    plot(t,Fce,'color',colors(i,:)); hold on
    plot(t,Fse,':','color',colors(i,:));
    plot(t,Fce+Fpe,'--','color',colors(i,:));
    
    subplot(235);
    plot(t,lce,'color',colors(i,:)); hold on
    
    subplot(236);
    plot(t,W,'color',colors(i,:)); hold on
    
    % summary measures
    Lmin(i) = min(lce);
    Lmax(i) = max(lce);
    Fce_max(i) = max(Fce);
    Fse_max(i) = max(Fse);

end

% make nice
titles = {'Excitation', 'Calcium activation', 'Force activation', 'Force','Length','Work'};
ylabels = {'u (0-1)', 'a_{Ca} (0-1)', 'a_F (0-1)', 'Force (N)', 'Length (m)', 'Work (J)'};

for i = 1:6
    subplot(2,3,i);
    box off
    title(titles{i})
    xlabel('Time (s)')
    ylabel(ylabels{i})
end

set(gcf,'units','normalized','position', [0 .3 .6 .4])

% summary figure
if ishandle(9), close(9); end; figure(9)
subplot(221)
plot(phis, Lmin,'linewidth',2); hold on
plot(phis, Lmax,'linewidth',2)

subplot(222);
plot(phis, Lmax-Lmin,'linewidth',2); hold on

subplot(223);
plot(phis, Fce_max*parms.mtc.r,'linewidth',2); hold on
plot(phis, Fse_max*parms.mtc.r,'linewidth',2)

subplot(224);
plot(Lmin, Fce_max,'linewidth',2); hold on
plot(Lmin, Fse_max,'linewidth',2)

% data from parms.mat
subplot(221)
plot(fl.phis, fl.Lces(:,1),'k--')
plot(fl.phis, fl.Lces(:,2),'k:')

subplot(222);
plot(fl.phis, -diff(fl.Lces,[],2),'k--')

subplot(223)
plot(fl.phis, fl.Fces(:,2)*parms.mtc.r,'k--')
plot(fl.phis, fl.Fses(:,2)*parms.mtc.r,'k:')

subplot(224)
plot(fl.Lces(:,2), fl.Fces(:,2),'k--')
plot(fl.Lces(:,2), fl.Fses(:,2),'k:');

titles = {'Length-angle','\DeltaLength-angle','Torque-angle','Force-length'}; %'Energy-angle','Energy-length'};
xlabels = {'Angle (deg)', 'Angle (deg)', 'Angle (deg)', 'Length (m)'};
ylabels = {'Length (m)', '\DeltaLength (m)', 'Torque (N-m)', 'Force (N)'};

for i = 1:4
    subplot(2,2,i); box off;
    title(titles{i})
    xlabel(xlabels{i})
    ylabel(ylabels{i})
    
    if i < 4
        xlim([0 max(phis)])
    end
end

set(gcf,'units','normalized','position', [.6 .3 .2 .4])

if saveFig
    save(fullfile(savedir,'Fig8','Fig8_force-length.mat'),'Lmin','Fce_max','Fse_max')
end

