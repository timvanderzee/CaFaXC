clear all; close all; clc

load('quad_parms.mat','parms','fl','fv');
parms.exp.phi = 90; % joint angle (full extension = 0)

parms.CB.f = 700;
parms.CB.g = [350 1300];
parms.CB.k = 5000;
parms.CB.b = 5000;

parms.CB.dLcrit = 2.3;
parms.CB.c = 1;

parms.CB.ps = 1;
parms.CB.w = .25;
parms.CB.rate_scale = 1;

% gaussian and the its integrals
parms.CB.gaussian.G =   @(x, c) c(1) * exp(-(x-c(2)).^2 / c(3));
parms.CB.gaussian.IG{1} = @(x, c) -0.5 * sqrt(pi*c(3))*c(1)*erf((c(2)-x)/sqrt(c(3)));
parms.CB.gaussian.IG{2} = @(x, c, G, IG1) c(2) * IG1(x, c) - 0.5 * c(3) * G(x, c);
parms.CB.gaussian.IG{3} = @(x, c, G, IG1) 0.5 * (2 * c(2)^2 + c(3)) * IG1(x, c) - 0.5 * c(3) * (c(2) + x) .* G(x, c);

parms.CB.K = 50;

parms.type = 'crossbridge_new';
parms.CB.Xmax = [.5 .5 .5];

% update
parms.ce.amin = 1e-3; % minimal excitation
parms = cfxc.calc_x0(parms); 

%% simulate isometric
parms.exp.A = 1;
parms.set.optimum = 1;

% simulate isometric contraction
parms.set.no_tendon = 1;
parms.exp.tstop = .05;

[t1,x1] = ode113(@cfxc.sim_muscle, [0 parms.exp.tstop], ones(1,6)*1e-3, parms.set.odeopt, parms);

for i = 1:size(x1,2)
    subplot(2,3,i)
    plot(t1, x1(:,i))
end

%%
dx = cfxc.sim_muscle(t1(end), x1(end,:),parms)