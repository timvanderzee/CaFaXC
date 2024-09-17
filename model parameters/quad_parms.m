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

%% cross-bridge model: cross-bridge cycling dynamics
% rate functions (only used for numerical approach)
parms.CB.f_func = @(parms) parms.CB.f(1) .* (parms.CB.xi>0 & parms.CB.xi<=1) .* parms.CB.xi;
parms.CB.g_func = @(parms) parms.CB.g(1) .* (parms.CB.xi >= 0) .* parms.CB.xi + ...
                           parms.CB.g(2) .* (parms.CB.xi < 0) + ...
                           parms.CB.g(3) .* (parms.CB.xi >= 1) .* (parms.CB.xi - 1);

% strain vector (only used for numerical approach)
parms.CB.xi = linspace(-5,5,1000);

% scale rates for submaximal activation
parms.CB.f = 150; % initial value, will be fit
parms.CB.g = [150 1000 150];  % initial value, will be fit

% original crossbridge model needs discretization
parms.CB.nbins = 20000;
parms.CB.xi0 = linspace(-50,50,parms.CB.nbins);

%% fit  crossbridge model rates on Hill-type force-velocity relation
% estimate steady-state isometric
parms.CB.Xmax = [1/2 1/4 1/6];
parms.CB.analytical = 1;

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
phis = (0:2.5:120); % [deg], joint angles
fl = cfxc.evaluate_force_length(phis, parms, show_figs);

%% activation and force facilitation dynamics
parms.ce.tau = [.002 .05]; % [s], forward and backward activation dynamics
parms.ce.tauR = [.06 .02]; % [s], forward and backward force facilitation dynamics

%% define some default settings
% experimental conditions
parms.type = 'CaFaXC';
parms.exp.stim_type = 'u_func';
parms.exp.vmtc = 0;
parms.exp.phi = 90; % [deg], joint angle

% settings
parms.set.optimum = 0;
parms.set.no_tendon = 0;
parms.set.odeopt = odeset('maxstep',1e-3);
parms.set.sim_mtc = 0; % simulate MTC
parms.set.fixed_velocity = 0;

%% saving
if save_parms
    cd(fullpath)
    save('quad_parms.mat','parms','fv','fl')
    disp('Saved: quad_parms.mat')
end




