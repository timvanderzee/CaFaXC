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

[fv.FCB(:,1),fv.n,fv.FCB(:,2)] = cfxc.CB_force_velocity(fv.vHill, parms);


%% fit  crossbridge model rates on Hill-type force-velocity relation
% rates used in the paper
parms.CB.g = [140 1388 78];
parms.CB.f = 140;

% option to refit
fit_CB = [];
parms.type = 'crossbridge';
parms.CB.ratefunc_type = 'Zahalak1981';
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
parms.ce.tau = [.002 .05]; % [s], forward and backward activation dynamics
parms.ce.tauR = [.06 .02]; % [s], forward and backward force facilitation dynamics

% fit time constants on twitch and tetanus data (using Hill-type)
parms.type = 'CaFaXC';


parms.exp.phi = 90; % joint angle (full extension = 0)
parms.ce.amin = 1e-3; % minimal excitation
parms = cfxc.calc_x0(parms); 

%% saving
if save_parms
    cd(fullpath)
    save('quad_parms.mat','parms','fv','fl')
    disp('Saved: quad_parms.mat')
end
