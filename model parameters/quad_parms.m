clear all; close all; clc

[fullpath, name, ext] = fileparts(which('quad_parms.m')); % folder contains the present m-file
cd(fullpath)
cd ..
addpath(genpath(cd))

save_parms = 1;
show_figs = 1;

%% start with some general parameters and (anonymous) functions
parms = cfxc.gen_funcs();
parms = cfxc.gen_parms(parms);

%% Hill-type force-velocity
parms.ce.vmaxrel    = 6; % [lopt/s] maximal contraction velocity
parms.ce.Arel       = .2; % [] curvature parameter
parms.ce.Fasymp     = 1.5; % eccentric asymptote

% Hill-type force-velocity
FHill = linspace(0, 0.99 * parms.ce.Fasymp); % force vector
vHill = parms.func.fv(1,  FHill, 1, parms); % corresponding velocities

% resample
fv.vHill = linspace(-parms.ce.vmaxrel, parms.ce.vmaxrel/2);
fv.FHill = interp1(vHill, FHill, fv.vHill);

%% fit  crossbridge model rates on Hill-type force-velocity relation
% rates used in the paper
parms.CB.g = [140 1388 78];
parms.CB.f = 140;

% option to refit
fit_CB = [];
[parms, fv] = cfxc.fit_CB_on_Hill(parms, fv,fit_CB);
cfxc.compare_fv(fv, parms)

%% CE force-length
parms.ce.Fmax = 4779; % [N]
parms.ce.lceopt = 0.09; % [m]
parms.ce.thickness = 0.04; % [m], if you're doing pennation    
parms.CB.delta = parms.CB.Xmax(2) / parms.ce.Fmax;

%% series-elastic element
parms.see.lse0 = .205; % [m], SE slack length
parms.see.sesl = 31; % [], SE linear stiffness
parms.see.sexm = .1; % [], SE toe-region length
parms.see.sefm = .5; % [], SE toe-region force

% calcualte shape so that things are continous
shape_fun = @(sh,parms) (sh * parms.see.sefm - parms.see.sesl*parms.see.sexm) * exp(sh) + parms.see.sesl*parms.see.sexm;
parms.see.sesh = fzero(@(sh) shape_fun(sh,parms), 4);
                        
%% muscle-tendon complex
parms.mtc.r = .042; % [m], moment arm 
parms.ce.lce0 = .057; % [m], CE length at phi = 0, determines position torque-angle relation
parms.mtc.lmtc0 = parms.see.lse0 + parms.ce.lce0;  % MTC length at phi = 0

%% evaluate torque-angle relation
fl = cfxc.evaluate_force_length((0:2.5:120), parms, show_figs);

%% define some default settings
% default model type
parms.type = 'CaFaXC';

% experimental conditions
parms.exp.stim_type = 'u_func';
parms.exp.vmtc = 0;
parms.exp.phi = 90; % [deg], joint angle

% settings
parms.set.optimum = 0; % 1 if ignorning CE force-length
parms.set.no_tendon = 0; % 1 if ignorning tendon
parms.set.sim_mtc = 0; % 1 if lMTC is a state
parms.set.fixed_velocity = 0; % 

% solver settings
parms.set.odeopt = odeset('maxstep',1e-3);

%% activation and force facilitation dynamics
clc
% load('quad_parms.mat')
parms.ce.tau = [.002 .05]; % [s], forward and backward activation dynamics
parms.ce.tauR = [.06 .02]; % [s], forward and backward force facilitation dynamics

% experimental observations
TTP = .1; % [s] time-to-peak (twitch)
T2T = 0.21; % twitch-to-tetanus ratio
HRT = .08; % [s] half relaxation time (twitch)
TRT = .19; % [s] 90% rise time of a tetanus

% fit time constants on twitch and tetanus data (using Hill-type)
parms.type = 'CaFaXC';

parms = cfxc.calc_x0(parms); 
parms.exp.tstop = .005;
parms.exp.A = 1;

X0 = [parms.ce.amin parms.exp.l0];
X0 = [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end)];

if ishandle(10), close(10); end
figure(10)

tstops = [.3 .005];
for i = 1:2
    parms.exp.tstop = tstops(i);

    [t,x] = ode23s(@cfxc.sim_muscle, [0 .3], X0, parms.set.odeopt, parms);

    lce = x(:,end);
    lse = parms.exp.lmtc - lce;
    Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;

    if i == 1
        Fmin = min(Fse);
        Fmax = max(Fse)-Fmin;
    end
        Frel = (Fse-Fmin)/Fmax;
        
        plot(t, (Fse-Fmin)/Fmax); hold on
    
    if i == 1
        x1 = t(find(Frel > .9, 1));
    else
        x2 = t(Frel == max(Frel));
        y2 = max(Frel);
        x3 = max(t(Frel > .5*y2));
    end
end

plot(TRT, .9,'o')
plot(TTP, T2T,'o')
plot(TTP+HRT, T2T*.5,'o');


plot(x1, .9, 'x')
plot(x2, y2, 'x')
plot(x3, y2/2, 'x')


data = [TTP T2T TTP + HRT];
cost = sum(([x2 y2 x3] - data).^2);

%%
% parms.x0 = X0;
% r0 = [parms.ce.tau(2) parms.ce.tauR];
% R = fminsearch(@(X) costfun(X, data, parms), r0);
% 
% %%
% parms.ce.tau(2) = R(1);
% parms.ce.tauR = R(2:3);

%% saving
if save_parms
    cd(fullpath)
    save('quad_parms.mat','parms','fv','fl')
    disp('Saved: quad_parms.mat')
end

%%
function[cost] = costfun(r, Y, parms)

tstops = [.3 .005];

parms.ce.tau(2) = r(1);
parms.ce.tauR = r(2:3);

ti = linspace(0,.3,100);

for i = 1:2
    parms.exp.tstop = tstops(i);

    [t,x] = ode23s(@cfxc.sim_muscle, [0 .3], parms.x0, parms.set.odeopt, parms);

    lce = x(:,end);
    lse = parms.exp.lmtc - lce;
    Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;

    if i == 1
        Fmin = min(Fse);
        Fmax = max(Fse)-Fmin;
    end
    
    Frel = (Fse-Fmin)/Fmax;
    
    if i == 1
        x1 = t(find(Frel > .9, 1));
    else
        x2 = t(Frel == max(Frel));
        y2 = max(Frel);
        x3 = max(t(Frel > .5*y2));
    end
    
    Fi(i,:) = interp1(t, Frel, ti);
    
end

X = [x2 y2 x3];

cost = sum((X - Y).^2);

figure(100)
plot(ti, Fi, '-', X(1), X(2), 'x', X(3), .5*X(2), 'x', ...
                  Y(1), Y(2), 'o', Y(3), .5*Y(2), 'o');  

drawnow
    
end




