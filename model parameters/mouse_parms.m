clear all; close all; clc
[fullpath, name, ext] = fileparts(which('quad_parms.m')); % folder contains the present m-file
cd(fullpath)
cd ..
addpath(genpath(cd))

load('quad_parms.mat')

save_parms = 1;
show_figs = 0;

%% Length and force scaling
% experiemntal observations Hollingworth
parms.ce.Fmax = 0.0017; % [N]
lopt = 0.014; % [m]

% assumed force-velocity
parms.ce.vmaxrel = 10;
parms.ce.Arel = .3;

% length scaling
lscale = lopt / parms.ce.lceopt;
parms.ce.lceopt    = parms.ce.lceopt  * lscale;
parms.see.lse0      = parms.see.lse0    * lscale;
parms.mtc.r         = parms.mtc.r       * lscale;
parms.mtc.lmtc0     = parms.mtc.lmtc0   * lscale;

%% Evaluate torque-angle relation
fl = cfxc.evaluate_force_length(0:2.5:120, parms, show_figs);

%% Activation and facilitation dynamics
parms.ce.tau = [.002 .012]; % [s], forward and backward activation dynamics
parms.ce.tauR = [.035 .06]; % [s], forward and backward force facilitation dynamics

%% crossbridge parameters
parms.CB.s = 2.4e-6; % [m], sarcomere length
parms.CB.mu = 1; % assume uniform

%% fit remaining crossbridge model parameters on Hill-type force-velocity relation
% analytical expressions Hill-type and Huxley-type force-velocity relations
[parms, fv] = cfxc.fit_CB_on_Hill(parms);

%% for DM model: evaluate isometric behavior (slightly different from Huxley), and update parms
parms = cfxc.get_DM_Xmax(parms, show_figs);

%% for DM model: evaluate force-velocity behavior (slightly different from Huxley)
fv = cfxc.evaluate_DM(parms, fv, show_figs);
 
%% save
parms.type = 'facilitation';

if save_parms
    cd(fullpath)
    save('mouse_parms.mat')
    disp('Saved: mouse_parms.mat')
end

