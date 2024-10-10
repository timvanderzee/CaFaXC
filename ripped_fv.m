clear all; close all; clc
load('quad_parms.mat')

% parms.CB.f = 700;
% parms.CB.g = [350 1300];
parms.CB.k = 5000;
parms.CB.b = 5000;

parms.CB.dLcrit = 2.3;
parms.CB.c = 1;

parms.CB.ps = 1;
parms.CB.w = .25;
parms.CB.rate_scale = 1;

parms.CB.f_func = @(parms) parms.f * (parms.xi > (parms.ps-parms.w) & parms.xi < (parms.ps+parms.w));
parms.CB.g_func = @(parms) parms.g(1) + parms.g(2) * (parms.xi < 0);

parms.CB.k_func = @(parms) parms.k * (parms.xi > parms.dLcrit);
parms.CB.b_func = @(parms) parms.b * (parms.xi > (parms.ps-parms.w) & parms.xi < (parms.ps+parms.w));

% gaussian and the its integrals
parms.CB.gaussian.G =   @(x, c) c(1) * exp(-(x-c(2)).^2 / c(3));
parms.CB.gaussian.IG{1} = @(x, c) -0.5 * sqrt(pi*c(3))*c(1)*erf((c(2)-x)/sqrt(c(3)));
parms.CB.gaussian.IG{2} = @(x, c, G, IG1) c(2) * IG1(x, c) - 0.5 * c(3) * G(x, c);
parms.CB.gaussian.IG{3} = @(x, c, G, IG1) 0.5 * (2 * c(2)^2 + c(3)) * IG1(x, c) - 0.5 * c(3) * (c(2) + x) .* G(x, c);
parms.CB.gaussian.IG{4} = @(x, c, G, IG1) 0.5 * c(2) * (2 * c(2)^2 + 3*c(3)) * IG1(x, c) - 0.5 * c(3) * (c(2)^2 + c(2)*x + c(3) + x.^2) .* G(x, c);

parms.CB.K = 50;


%% determine force-velocity
if ishandle(10), close(10); end

parms.set.fixed_velocity = 1;

parms.type = 'crossbridge';

parms.CB.f = 454;
parms.CB.g = [454 2100 88];
parms.CB.Xmax = [1/2 1/4 1/6];

fv.vHill = linspace(-12,6,100);
[fv.FCB(:,1),fv.n,fv.FCB(:,2)] = cfxc.CB_force_velocity(fv.vHill, parms);

parms.CB.ratefunc_type = 'Zahalak1981';
[fv.FCB(:,3), parms] = cfxc.evaluate_DM(fv.vHill, parms);

parms.CB.f = 454;
parms.CB.g = [454 2100/2 88];
parms.CB.Xmax = [1/2 1/2 .6];

parms.CB.ratefunc_type = 'vanderZee2024';
[fv.FCB(:,4), parms] = cfxc.evaluate_DM(fv.vHill, parms);

parms.CB.dLcrit = 2.5;
parms.CB.Xmax = [1/2 1/2 .6 0];

parms.CB.k1 = 1;
parms.CB.kF = 50;
parms.CB.k2 = 10;

parms.CB.f = 1000;
parms.CB.g = [1000 2100/2 88];

parms.type = 'crossbridge_new';
[fv.FCB(:,5), parms] = cfxc.evaluate_DM(fv.vHill, parms);

figure(10)
plot(fv.vHill, fv.FCB)

legend('Huxley','Huxley v2', 'Zahalak', 'vanderZee', 'Ripped')

%%
A = 1;

vs = linspace(-10,5,10);

for j = 1:length(vs)
    disp(j)    
    u_imposed = vs(j) * 0.5*parms.CB.s / (parms.CB.h);

    % update
    parms.ce.amin = 1e-3; % minimal excitation
    parms.CB.Xmax = [.5 .5 .6];
    parms = cfxc.calc_x0(parms); 

    %% simulate isometric
    parms.exp.A = A;
    parms.set.optimum = 1;
    parms.set.fixed_velocity = 1;

    % simulate isometric contraction
    parms.set.no_tendon = 1;
    parms.exp.tstop = .05;
    
    X0 = ones(1,5)*1e-3;
    [t0,x0] = ode113(@cfxc.sim_muscle, [0 parms.exp.tstop], X0, parms.set.odeopt, parms);

    %% simulate stretch
    parms.exp.u = u_imposed;

    [t1,x1] = ode113(@cfxc.sim_muscle, [0 parms.exp.tstop], x0(end,:), parms.set.odeopt, parms);
    t = [t0; t1+t0(end)];
    x = [x0; x1];
% 
%     figure(1)
%     for i = 1:3
%         subplot(2,3,i)
%         plot(t, x(:,i+1)./parms.CB.Xmax(i)); hold on
%     end
%     
    Fss(j) = x(end,3)/parms.CB.Xmax(2);
end

%%
if ishandle(1), close(1); end; figure(1)
plot(fv.vHill, fv.FCB); hold on

figure(1)
plot(vs, Fss,'o');


