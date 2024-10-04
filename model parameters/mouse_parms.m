clear all; close all; clc

[fullpath, name, ext] = fileparts(which('quad_parms.m')); % folder contains the present m-file
cd(fullpath)
cd ..
addpath(genpath(cd))

save_parms = 1;
show_figs = 1;

%% Length and force scaling
load('quad_parms.mat') % start with quad parms

% experimental observations of Hollingworth
parms.ce.Fmax   = .0017; % [N]
lopt            = .014; % [m]

% assumed force-velocity
parms.ce.vmaxrel    = 10;
parms.ce.Arel       = .3;

% Hill-type force-velocity
FHill = linspace(0, 0.99 * parms.ce.Fasymp); % force vector
vHill = parms.func.fv(1,  FHill, 1, parms); % corresponding velocities

% resample
fv.vHill = linspace(-parms.ce.vmaxrel, parms.ce.vmaxrel/2);
fv.FHill = interp1(vHill, FHill, fv.vHill);

% length scaling
lscale = lopt / parms.ce.lceopt;
parms.ce.lceopt    = parms.ce.lceopt  * lscale;
parms.see.lse0      = parms.see.lse0    * lscale;
parms.mtc.r         = parms.mtc.r       * lscale;
parms.mtc.lmtc0     = parms.mtc.lmtc0   * lscale;

%% crossbridge parameters
parms.CB.s = 2.4e-6; % [m], sarcomere length
parms.CB.mu = 1; % assume uniform fiber type

%% fit  crossbridge model rates on Hill-type force-velocity relation
% rates used in the paper
parms.CB.g = [287 2049 161];
parms.CB.f = 287;

% option to refit
[parms, fv] = cfxc.fit_CB_on_Hill(parms, fv,'fit_CB');
cfxc.compare_fv(fv, parms)

%% Evaluate torque-angle relation
fl = cfxc.evaluate_force_length(0:2.5:120, parms, show_figs);

%% Activation and facilitation dynamics
parms.ce.tau = [.002 .012]; % [s], forward and backward activation dynamics
parms.ce.tauR = [.035 .06]; % [s], forward and backward force facilitation dynamics
 
%% save
if save_parms
    cd(fullpath)
    save('mouse_parms.mat')
    disp('Saved: mouse_parms.mat')
end

