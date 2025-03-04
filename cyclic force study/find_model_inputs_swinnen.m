clear all; close all; clc

save_output = 1;

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
Tamp = 30 - Tmin; % [N-m]

%% load parameters
addpath(fullfile(cfxcpath,'model parameters'));

muscle = 'calf';
load([muscle,'_parms.mat'])

% simulate minimal and maximal
parms.type = 'crossbridge';

parms.exp.phi = 0;
parms = cfxc.calc_x0(parms);

%% define kinematics and kinetics
freqs = 0.5:0.1:4;
freqs = .5:.5:2.5;
% freqs = 0.5;

amps = [0 5 10];
freqs = 2 * ones(size(amps));

phis = 20 * ones(length(freqs), N);
tx = zeros(length(freqs),N);

for f = length(freqs):-1:1
    freq = freqs(f);
    
    tx(f,:) = linspace(0, 1/freq, N);
    Tx = Tamp/2 - Tamp/2 * cos(2*pi*freq*tx(f,:)) + Tmin;
    
    phis(f,:) = amps(f) * cos(freqs(f)*2*pi*tx(f,:)) - 10;

end

lmtcs = parms.func.lmtc(phis, parms);

figure(1)
subplot(121); plot(tx', Tx'); title('Kinetics (torque)'); ylabel('Torque (N-m)')
subplot(122); plot(tx', phis'); title('Kinematics (angle)'); ylabel('Angle (deg)')

for i = 1:2
    subplot(1,2,i); 
    xlabel('Time (s)'); 
    box off; 
%     ylim([0 25]);
end

%% find muscle length and velocity
Fx = Tx ./ parms.mtc.r;
lses = parms.see.lse0 * parms.func.lse(Fx(f,:)/parms.ce.Fmax, parms) + parms.see.lse0;
lces = lmtcs - lses;
FL = parms.func.fce(lces, parms);
    
vce = nan(length(freqs), N);
for f = length(freqs):-1:1
    vce(f,:) = cfxc.grad5(lces(f,:)', mean(diff(tx(f,:))))'; % positive for shortening
end

figure(2)
subplot(131); plot(tx', lces');
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
Ux = vce * parms.CB.s / (2*parms.CB.h*parms.ce.lceopt);

% create problem
prob.prev = 1;
prob.N = N;
prob.freqs = freqs;
prob.M = 3;
prob.Target = repmat(Fx(:)/parms.ce.Fmax, 1, length(freqs),1);
prob.t = tx;

if save_output
save(fullfile(SAVEDIR,'kinematics_kinetics_swinnen.mat'))
end

%% find cross-bridge activation
parms.Xmin = parms.exp.x0(2:4);
parms.Umin = parms.ce.amin;

% create parms matrix
Parms = parms;
for f = 1:length(freqs)
    for i = 1:N
        Parms(i,f) = parms(1);

        Parms(i,f).exp.u = Ux(f,i);

        % force-length
        Parms(i,f).exp.a = FL(f,i);
    end
end

[R_opt, X_opt] = cfxc.solve_optimal_control(prob, Parms, SAVEDIR);

% save?
if save_output
save_xr_name = fullfile(SAVEDIR,'crossbridge_activations_all_swinnen.mat');
save(save_xr_name,'X_opt','R_opt');
end

%% find facilitation 
parms.type = 'facilitation';
prob.prev = 0;
prob.M = 1;
prob.Target = R_opt';
prob.N = length(tx);
prob.freqs = freqs;


[~,ipeak] = max(R_opt,[],2);

% create parms matrix
% parms.tau = [];
Parms = parms;
for f = 1:length(freqs)

Rdot = cfxc.grad5(R_opt(f,:)', 1);
    
    for i = 1:prob.N
        Parms(i,f) = parms(1);
        if Rdot(i) > 0
            Parms(i,f).ce.tau = parms.ce.tauR(1);
        else
            Parms(i,f).ce.tau = parms.ce.tauR(2);
        end
    end
end

[C_opt, R_opt] = cfxc.solve_optimal_control(prob, Parms, SAVEDIR); 
[C_opt_smooth] = cfxc.smoothen_solution(C_opt, prob);

if save_output
    save_ca_name = fullfile(SAVEDIR,'calcium_activations_swinnen.mat');
    save(save_ca_name,'C_opt','R_opt','C_opt_smooth');
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
Cdot = cfxc.grad5(C_opt_smooth(f,:)', 1);

    for i = 1:prob.N
        Parms(i,f) = parms(1);
        if Cdot(i) > 0
            Parms(i,f).ce.tau = parms.ce.tau(1);
        else
            Parms(i,f).ce.tau = parms.ce.tau(2);
        end
    end
end

[U_opt, C_opt] = cfxc.solve_optimal_control(prob, Parms); 
[U_opt_smooth] = cfxc.smoothen_solution(U_opt, prob);

if save_output
    save_me_name = fullfile(SAVEDIR,'muscle_excitations_swinnen.mat');
    save(save_me_name,'U_opt','C_opt','U_opt_smooth')
end





