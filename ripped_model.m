clear all; close all; clc

load('quad_parms.mat','parms','fl','fv');
parms.exp.phi = 90; % joint angle (full extension = 0)

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
parms.CB.gaussian.IG{4} = @(x, c, G, IG1) 0.5 * c(2) * (2 * c(2)^2 + 3*c(3)) * IG1(x, c) - 0.5 * c(3) * (c(2)^2 + c(2)*x + c(3) + x.^2) .* G(x, c);

parms.CB.K = 50;
close all


%% activation
parms.type = 'crossbridge';
parms.CB.Xmax = [.5 .5 .5];
parms.CB.f = 350;
parms.CB.g = [350 1300 80];
X0 = ones(1,4)*1e-3;

parms.CB.ratefunc_type = 'vanderZee2024';

% parms.type = 'crossbridge';
% parms.CB.Xmax = [1/2 1/4 1/6];
% parms.CB.f = 454;
% parms.CB.g = [454 2100 88];
% X0 = ones(1,4)*1e-3;


parms.exp.A = .5;
parms.set.optimum = 1;
parms.set.fixed_velocity = 1;

% simulate isometric contraction
parms.set.no_tendon = 1;
parms.exp.tstop = .05;

% update
parms.ce.amin = 1e-3; % minimal excitation
parms = cfxc.calc_x0(parms); 

As = linspace(0,1,10);

for i = 1:length(As)
    parms.exp.A = As(i);
    [t,x] = ode113(@cfxc.sim_muscle, [0 parms.exp.tstop], X0, parms.set.odeopt, parms);

    figure(1)
    for j = 1:3
        subplot(1,3,j)
        plot(t, x(:,j+1)./parms.CB.Xmax(j)); hold on
    end
end

%%
close all
u_imposed = 500;
A = 1;

parms.CB.K = 10;

figure(1)
color = get(gca,'colororder');

for j = 1:4
    if j == 1
        
        parms.CB.f_func = @(parms) parms.f(1) .* (parms.xi>0 & parms.xi<=1) .* parms.xi;
    
        parms.CB.g_func = @(parms) parms.g(1) .* (parms.xi >= 0) .* parms.xi + ...
                           parms.g(2) .* (parms.xi < 0) + ...
                           parms.g(3) .* (parms.xi >= 1) .* (parms.xi - 1);
        parms.type = 'crossbridge';
        parms.CB.Xmax = [1/2 1/4 1/6];
        parms.CB.f = 454;
        parms.CB.g = [454 2100 88];
        X0 = ones(1,4)*1e-3;
        
        parms.CB.analytical = 1;
        parms.CB.ratefunc_type = 'Zahalak1981';
        
    elseif j == 2
        parms.CB.analytical = 0;
        
    elseif j == 3
        parms.type = 'crossbridge';
        parms.CB.f = 454;
        parms.CB.g = [454 2100/2 88];
        parms.CB.Xmax = [.5 .5 .6];

        X0 = ones(1,4)*1e-3;
        parms.CB.f_func = @(parms) parms.f * (parms.xi > (parms.ps-parms.w) & parms.xi < (parms.ps+parms.w));
        parms.CB.g_func = @(parms) parms.g(1) * (parms.xi > 0 & parms.xi < 1) + parms.g(2) * (parms.xi < 0) + ((parms.g(3)).*parms.xi + parms.g(1) - parms.g(3)) .* (parms.xi > 1);

        parms.CB.analytical = 1;
        parms.CB.ratefunc_type = 'vanderZee2024';
        
    elseif j == 4
        parms.type = 'crossbridge_new';
%         parms.CB.Xmax = [.5 .5 .5];
%         parms.CB.f = 700;
%         parms.CB.g = [350 1300];
        X0 = ones(1,5)*1e-3;
        parms.CB.analytical = 1;
        
    elseif j == 5
        parms.type = 'crossbridge_new';
%         parms.CB.Xmax = [.5 .5 .5];
%         parms.CB.f = 700;
%         parms.CB.g = [350 1300];
        X0 = ones(1,5)*1e-3;
        parms.CB.analytical = 1;
    end
    
%     X0(4) = X0(3)*1.5;
    % update
    parms.ce.amin = 1e-3; % minimal excitation
    parms = cfxc.calc_x0(parms); 

    %% simulate isometric
    parms.exp.A = A;
    parms.set.optimum = 1;
    parms.set.fixed_velocity = 1;

    % simulate isometric contraction
    parms.set.no_tendon = 1;
    parms.exp.tstop = .05;

    [t0,x0] = ode113(@cfxc.sim_muscle, [0 parms.exp.tstop], X0, parms.set.odeopt, parms);

    %% simulate stretch
    parms.exp.u = u_imposed;

    [t1,x1] = ode113(@cfxc.sim_muscle, [0 parms.exp.tstop], x0(end,:), parms.set.odeopt, parms);
    t = [t0; t1+t0(end)];
    x = [x0; x1];

    figure(1)
    for i = 1:3
        subplot(2,3,i)
        plot(t, x(:,i+1)./parms.CB.Xmax(i)); hold on
    end
    Fss = x(end,3)/parms.CB.Xmax(2);
    
    for i = 1:3
        subplot(2,3,i+3)
        plot(t, x(:,i+1)); hold on
    end
    
    n0 = cfxc.n_func(x0(end,2), x0(end,3), x0(end,4), parms.CB);
    n1 = cfxc.n_func(x1(end,2), x1(end,3), x1(end,4), parms.CB);
    
    figure(2)
    plot(parms.CB.xi, n0,'--','color',color(j,:)); hold on
    plot(parms.CB.xi, n1,'color',color(j,:)); hold on
    
end

%%
figure(1)
legend('Zahalak ana','Zahalak num','vanderZee ana','vanderZee num','ripping')

figure(2)
legend('Zahalak ana','-','Zahalak num', '-','vanderZee ana','-','vanderZee num', '-','ripping','-')

figure(3)
parms.type = 'crossbridge_new';


fv.vHill = linspace(-12,6,100);
[fv.FCB, parms] = cfxc.evaluate_DM(fv.vHill, parms);


figure(3)
plot(fv.vHill, fv.FCB, '-', u_imposed*parms.CB.h/(.5*parms.CB.s), Fss,'o')


%% compare dx 
parms.type = 'crossbridge';
X0 = ones(1,4)*1e-3;
parms.CB.analytical = 1;
parms.CB.ratefunc_type = 'vanderZee2024';


dx = cfxc.sim_muscle(0, [1 .5 .5 .6], parms)





