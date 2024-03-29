clear all; close all; clc

save_output = 0;

%% constants for dependencies
[cfxcpath, name, ext] = fileparts(which('main_cfxc.m')); % folder contains the present m-file
addpath(genpath(cfxcpath));

import casadi.*

%% set up constants. 

% save directoy
SAVEDIR = fullfile(cfxcpath,'cyclic force study','inputs');
% nodes for collocation.
N = 100;

% target properties
Tmin = .2; % [N-m], mimimal torque
Tamp = 11 - Tmin; % [N-m]

%% load parameters
addpath(fullfile(cfxcpath,'model parameters'));

muscle = 'quad';
load([muscle,'_parms.mat'])

% simulate minimal and maximal
parms.type = 'crossbridge';
parms = cfxc.update_parms(parms);

%% define kinematics and kinetics
freqs = 0.5:0.1:2.5;
phi = 20 * ones(length(freqs), N);
lmtc = parms.lmtc_func(phi, parms);

for f = length(freqs):-1:1
    freq = freqs(f);
    
    tx(f,:) = linspace(0, 1/freq, N);
    Tx = Tamp/2 - Tamp/2 * cos(2*pi*freq*tx(f,:)) + Tmin;
end

figure(1)
subplot(121); plot(tx', Tx'); title('Kinetics (torque)'); ylabel('Torque (N-m)')
subplot(122); plot(tx, phi); title('Kinematics (angle)'); ylabel('Angle (deg)')

for i = 1:2
    subplot(1,2,i); 
    xlabel('Time (s)'); 
    box off; 
    ylim([0 25]);
end

%% find muscle length and velocity
Fx = Tx ./ parms.r;
lse = parms.lse0 * parms.lse_func(Fx(f,:)/parms.Fmax, parms) + parms.lse0;
lce = lmtc - lse;
FL = parms.afunc(lce, parms);
    
for f = length(freqs):-1:1
    vce(f,:) = cfxc.grad5(lce(f,:)', mean(diff(tx(f,:))))'; % positive for shortening
end

figure(2)
subplot(131); plot(tx', lce');
subplot(132); plot(tx', FL');
subplot(133); plot(tx', vce');

titles = {'l_{CE}', 'f_L','v_{CE}'};
units = {'m','','m/s'};
for i = 1:3
    subplot(1,3,i); 
    xlabel('Time (s)'); 
    box off; 
    title(titles{i})
    ylabel([titles{i}, ' (',units{i},')']);
end

% scale velocity
Ux = vce * parms.s / (2*parms.h*parms.lceopt);

% create problem
prob.prev = 1;
prob.N = N;
prob.freqs = freqs;
prob.M = 3;
prob.Target = repmat(Fx(:)/parms.Fmax, 1, length(freqs),1);
prob.t = tx;

if save_output
save(fullfile(SAVEDIR,'kinematics_kinetics.mat'))
end

%% find cross-bridge activation
parms.Xmin = parms.x0(2:4);
parms.Umin = parms.amin;

% create parms matrix
Parms = parms;
for f = 1:length(freqs)
    for i = 1:N
        Parms(i,f) = parms(1);

        Parms(i,f).u = Ux(f,i);

        % force-length
        Parms(i,f).a = FL(f,i);
    end
end

[R_opt, X_opt] = cfxc.solve_optimal_control(prob, Parms, SAVEDIR);

% save?
if save_output
save_xr_name = fullfile(SAVEDIR,'crossbridge_activations_all.mat');
save(save_xr_name,'X_opt','R_opt');
end
% load(['crossbridge_activations.mat'],'X_opt','R_opt')
% R_opt = Rc;

%% find facilitation 
parms.type = 'facilitation';
prob.prev = 1;
prob.M = 1;
prob.Target = R_opt';
prob.N = length(tx);
prob.freqs = freqs;

[~,ipeak] = max(R_opt,[],2);

% create parms matrix
Parms = parms;
for f = 1:length(freqs)
    for i = 1:prob.N
        Parms(i,f) = parms(1);
        if i < ipeak(f) 
            Parms(i,f).tau = parms.tauR(1);
        else
            Parms(i,f).tau = parms.tauR(2);
        end
    end
end

[C_opt, R_opt] = cfxc.solve_optimal_control(prob, Parms, SAVEDIR); 
[C_opt_smooth] = cfxc.smoothen_solution(C_opt, prob);

if save_output
save('calcium_activations.mat','C_opt','R_opt','C_opt_smooth');
end

%% find excitation 
% load('calcium_activations.mat','R_opt','C_opt_smooth');
parms.type = 'activation';
prob.prev = 0;
prob.M = 1;
prob.Target = C_opt_smooth';

% parms.tau(1) = 0.01;

[~,ipeak] = max(C_opt_smooth,[],2);

% create parms matrix
Parms = parms;
for f = 1:length(freqs)
    for i = 1:prob.N
        Parms(i,f) = parms(1);
        if i < ipeak(f) 
            Parms(i,f).tau = parms.tau(1);
        else
            Parms(i,f).tau = parms.tau(2);
        end
    end
end

[U_opt, C_opt] = cfxc.solve_optimal_control(prob, Parms); 
[U_opt_smooth] = cfxc.smoothen_solution(U_opt, prob);

if save_output
save('muscle_excitations.mat','U_opt','C_opt','U_opt_smooth')
end





