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
parms.ce.vmaxrel = 6; % [lopt/s] maximal contraction velocity
parms.ce.Arel = .2; % [] curvature parameter

% Hill-type force-velocity
FHill = linspace(0, 0.99 * parms.ce.Fasymp); % force vector
vHill = parms.func.fv(1,  FHill, 1, parms); % corresponding velocities

% resample
fv.vHill = linspace(-parms.ce.vmaxrel, parms.ce.vmaxrel/2);
fv.FHill = interp1(vHill, FHill, fv.vHill);

%% fit  crossbridge model rates on Hill-type force-velocity relation
% analytical expressions Hill-type and Huxley-type force-velocity relations
[parms, fv] = cfxc.fit_CB_on_Hill(parms, fv);

%% CE force-length
parms.ce.Fmax = 4779; % [N]
parms.ce.lceopt = 0.093; % [m]
parms.ce.thickness = 0.04; % [m], if you're doing pennation    
parms.CB.delta = parms.CB.Xmax(2) / parms.ce.Fmax;

%% series-elastic element
parms.see.lse0 = .205; % [m], SE slack length
parms.see.sexm = .097; % [], SE toe-region length
parms.see.sesl = 31; % [], SE linear stiffness
parms.see.sefm = 0.5; % [], SE toe-region force

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
parms.ce.tau = [.002 .05]; % [s], forward and backward activation dynamics
parms.ce.tauR = [.06 .02]; % [s], forward and backward force facilitation dynamics

% experimental observations
TTP = .1; % [s] time-to-peak (twitch)
T2T = 0.21; % twitch-to-tetanus ratio
HRT = .1; % [s] half relaxation time (twitch)

% fit time constants on twitch and tetanus data (using Hill-type)
parms.type = 'CaFaXC';

parms = cfxc.calc_x0(parms); 
parms.exp.tstop = .005;
parms.exp.A = 1;

X0 = [parms.ce.amin parms.exp.l0];
X0 = [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end)];
[t,x] = ode23s(@cfxc.sim_muscle, [0 1], X0, parms.set.odeopt, parms);
 
lce = x(:,end);
lse = parms.exp.lmtc - lce;
Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;

%%
if ishandle(10), close(10); end
figure(10)
plot(t, Fse/parms.ce.Fmax); hold on
plot(TTP, T2T,'o')
plot(TTP*2, T2T*.5,'o');

%% saving
if save_parms
    cd(fullpath)
    save('quad_parms_v2.mat','parms','fv','fl')
    disp('Saved: quad_parms.mat')
end




