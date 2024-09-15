clear all; close all; clc

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

%% Do a twitch
close all
clc
parms.exp.tstop = .005;
parms.exp.A = 1;

color = get(gca,'colororder');

for j = 1:3
    if j == 1
        parms.type = 'CaFaXC';
        X0 = [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end)];
    elseif j == 2
        parms.type = 'CaFaXC';
        X0 = [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end-1)];   
    else
        parms.type = 'CaFaXC_v2';
        X0 = [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2) parms.exp.x0(4) parms.exp.x0(end)];  
    end
    
    [t,x] = ode23s(@cfxc.sim_muscle, [0 1], X0, parms.set.odeopt, parms);

    X = x;
    if j == 3
        X(:,6) = x(:,5);
        X(:,5) = x(:,4);
        X(:,4) = nan;
    end
    
    figure(1)
    for i = 1:size(X,2)
        subplot(3,2,i)
        plot(t,X(:,i),'color',color(j,:)); hold on
    end
end
