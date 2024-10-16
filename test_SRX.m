clear all; close all; clc

load('quad_parms.mat');
parms.CB.k = 100;
parms.CB.b = 5000;

parms.CB.dLcrit = 2.3;
parms.CB.c = 1;

parms.CB.ps = 1;
parms.CB.w = .25;
parms.CB.rate_scale = 1;

parms.CB.g_func = @(parms) parms.g(1) + parms.g(2) * (parms.xi < 0);
parms.CB.k_func = @(parms) parms.k * (parms.xi > parms.dLcrit);

parms.CB.b_func = @(parms) parms.b * (parms.xi > (parms.ps-parms.w) & parms.xi < (parms.ps+parms.w));
parms.CB.f_func = @(parms) parms.f * (parms.xi > (parms.ps-parms.w) & parms.xi < (parms.ps+parms.w));

% gaussian and the its integrals
parms.CB.gaussian.G =   @(x, c) c(1) * exp(-(x-c(2)).^2 / c(3));
parms.CB.gaussian.IG{1} = @(x, c) -0.5 * sqrt(pi*c(3))*c(1)*erf((c(2)-x)/sqrt(c(3)));
parms.CB.gaussian.IG{2} = @(x, c, G, IG1) c(2) * IG1(x, c) - 0.5 * c(3) * G(x, c);
parms.CB.gaussian.IG{3} = @(x, c, G, IG1) 0.5 * (2 * c(2)^2 + c(3)) * IG1(x, c) - 0.5 * c(3) * (c(2) + x) .* G(x, c);
parms.CB.gaussian.IG{4} = @(x, c, G, IG1) 0.5 * c(2) * (2 * c(2)^2 + 3*c(3)) * IG1(x, c) - 0.5 * c(3) * (c(2)^2 + c(2)*x + c(3) + x.^2) .* G(x, c);

parms.CB.K = 200;
% parms.CB.k1 = .5;
% parms.CB.kF = 25;
% parms.CB.k2 = 5;



%% isometric
parms.thin_filament_coop = 1;
parms.CB.K = 1e5;
parms.set.optimum = 1;

clc
parms.type = 'crossbridge_new';
% parms.CB.ratefunc_type = 'vanderZee2024';
parms.set.fixed_velocity = 1;
parms.exp.u = 0;

parms.CB.kcoop = 50;
parms.CB.kon = 10;
parms.CB.koff = 10;

parms.exp.x0 = [0 .01 .01 .015 0 0];
parms.set.sim_mtc = 0;
parms.exp.A = 1;
parms.exp.tstop = 3;

parms.CB.f = 140;
[t0,x0] = ode113(@cfxc.sim_muscle, [0 3], parms.exp.x0, parms.set.odeopt, parms);

% reset XB distribution
parms.CB.f = 0;
X0 = parms.exp.x0;
X0(end) = x0(end,end);
tic
[t1,x1] = ode113(@cfxc.sim_muscle, [0 3], X0, parms.set.odeopt, parms);
toc
% reset XB distribution
parms.CB.f = 140;
[t2,x2] = ode113(@cfxc.sim_muscle, [0 3], x1(end,:), parms.set.odeopt, parms);

t = [t0; t1+t0(end); t2+t1(end)+t0(end)];
x = [x0; x1; x2];

close all
figure(1)
for i = 1:size(x,2)
    subplot(2,3,i)
    plot(t, x(:,i)); hold on
end

%% simulate muscle stretch
parms.type = 'crossbridge_new';
parms.CB.k = 2000;
parms.CB.b = 5000;
parms.set.optimum = 1;

% parms.CB.Xmax = [.5 .5 .51 0 1];
% parms.CB.Xmax = [.5 .5 .51];
% parms = cfxc.calc_x0(parms); 

parms.exp.tstop = 5;
parms.set.fixed_velocity = 1;

v = 500;
vs = [0 v -v v];
ts = [5 1 2 1];
as = [1 1 1 1];

% parms.exp.x0 = parms.exp.x0(1:end-1);

[t,x] = cfxc.stretch_protocol(as, vs, ts, parms, 'muscle');

close all
figure(1)
for i = 1:size(x,2)
    subplot(2,3,i)
    plot(t, x(:,i)); hold on
    box off
end

%% determine initial conditions
parms.type = 'crossbridge_new';
parms.exp.tstop = 5;
parms.set.fixed_velocity = 0;
parms.exp.phi = 70;

parms.CB.Xmax = [.5 .5 .51 0 .5];
% parms.CB.Xmax = [.5 .5 .51];

parms.exp.phi = 0;
parms = cfxc.calc_x0(parms);     
X0 = parms.exp.x0;

%% simulate muscle-tendon stretch
% parms.exp.tstop = 
vs = [0 1 -1 1] * .1;
ts = [2 1 1 1];

parms.CB.kon = 5;
parms.CB.koff = 5;

parms.set.optimum = 0;
parms.set.sim_mtc = 1;

parms.exp.x0 = [X0 parms.exp.lmtc];
    
As = [.5 1];
close all    

for j = 1:length(As)
    as = As(j) * ones(size(vs));

    [t,x] = cfxc.stretch_protocol(as, vs, ts, parms, 'muscle-tendon');
    [y,X] = cfxc.get_model_output(t, x, parms);

   
    figure(1)
    titles = {'Calcium','Q0','Q1','Q2','R','N','L_M','L_{MT}'};

    for i = 1:size(x,2)
        subplot(2,4,i)
        plot(t, x(:,i)); hold on

        box off
    %     xlim([ts(1) max(t)])
        title(titles{i})
    end
end

%% initial conditions
parms.exp.phi = 20;
phi_rad_0 = -parms.exp.phi * (pi/180);

parms.exp.stim_type = 'constant';
parms.exp.A = .1;

parms.type = 'crossbridge_new';
parms.set.optimum = 0; %
parms.set.sim_mtc = 0; %

% parms has too many parameters; here we specify explicity for some code control.
parms.CB.Xmax = [.5 .5 .51 0 .5];
parms = cfxc.calc_x0(parms); 

% simule max contraction
parms.set.optimum = 0;
parms.exp.A = .01;
[t,x] = ode113(@cfxc.sim_muscle, [0 5], parms.exp.x0, parms.set.odeopt, parms);

close all
figure(1)
plot(t,x)

X0 = x(end,:);

%% sim leg
parms.set.fixed_velocity = 0;
parms.set.optimum = 0; %
parms.set.sim_mtc = 0; % 

close all
parms.mode = 'pre-movement';

% guess from Willaert et al. (2024) Fig. 1
dphi = 40; % [deg]
parms.f = 0.5;
T = 1/parms.f;

parms.exp.A = .01;

parms.CB.f = 20;
parms.CB.g = [20 100];

for j = 1:2
     X0f = [X0, -parms.exp.phi, 0];
%     if j == 1
%         % add angle and angular velocity
%         X0f = [parms.exp.x0, phi_rad_0, 0];
%     else
%         X0f = [X0, phi_rad_0, 0];
%     end

    if j == 1
        parms.A = dphi/2;
    else
        parms.A = .001;
    end

    parms.mode = 'pre-movement';
    [t0,x0] = ode113(@cfxc.sim_segment, [0 3*T], X0f, [], parms);

    parms.mode = 'regular';
    [t1,x1] = ode113(@cfxc.sim_segment, [0 5], x0(end,:), [], parms);
    
    t = [t0; t1+t0(end)];
    x = [x0; x1];
    
    figure(1)
    titles = {'Calcium','Q0','Q1','Q2','R','N','L_M','phi','omega'};

    for i = 1:size(x,2)
        subplot(2,5,i)
        plot(t, x(:,i)); hold on

        box off
        title(titles{i})
        xline(3*T,'k--')

    end
end

%%
dx = cfxc.sim_segment(0,x0(end,:)', parms);
% clcs = cfxc.calcs_facilitation_sim(t_s,state_s,parms);
