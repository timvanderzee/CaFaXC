clc; close all; clear all; 

% experimental observations
TTP = .1; % [s] time-to-peak (twitch)
T2T = 0.21; % twitch-to-tetanus ratio
HRT = .08; % [s] half relaxation time (twitch)
TRT = .19; % [s] 90% rise time of a tetanus

load('quad_parms.mat')

parms = cfxc.calc_x0(parms); 
parms.exp.tstop = .005;
parms.exp.A = 1;

X0 = [parms.ce.amin parms.exp.l0];
X0 = [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end)];

if ishandle(10), close(10); end
figure(10)

tstops = [.3 .005];
for i = 1:2
    parms.exp.tstop = tstops(i);

    [t,x] = ode23s(@cfxc.sim_muscle, [0 .3], X0, parms.set.odeopt, parms);

    lce = x(:,end);
    lse = parms.exp.lmtc - lce;
    Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;

    if i == 1
        Fmin = min(Fse);
        Fmax = max(Fse)-Fmin;
    end
        Frel = (Fse-Fmin)/Fmax;
        
        plot(t, (Fse-Fmin)/Fmax); hold on
    
    if i == 1
        x1 = t(find(Frel > .9, 1));
    else
        x2 = t(Frel == max(Frel));
        y2 = max(Frel);
        x3 = max(t(Frel > .5*y2));
    end
end

plot(TRT, .9,'o')
plot(TTP, T2T,'o')
plot(TTP+HRT, T2T*.5,'o');


plot(x1, .9, 'x')
plot(x2, y2, 'x')
plot(x3, y2/2, 'x')


data = [TTP T2T TTP + HRT];
cost = sum(([x2 y2 x3] - data).^2);

%%
% parms.x0 = X0;
% r0 = [parms.ce.tau(2) parms.ce.tauR];
% R = fminsearch(@(X) costfun(X, data, parms), r0);
% 
% %%
% parms.ce.tau(2) = R(1);
% parms.ce.tauR = R(2:3);

function[cost] = costfun(r, Y, parms)

tstops = [.3 .005];

parms.ce.tau(2) = r(1);
parms.ce.tauR = r(2:3);

ti = linspace(0,.3,100);

for i = 1:2
    parms.exp.tstop = tstops(i);

    [t,x] = ode23s(@cfxc.sim_muscle, [0 .3], parms.x0, parms.set.odeopt, parms);

    lce = x(:,end);
    lse = parms.exp.lmtc - lce;
    Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;

    if i == 1
        Fmin = min(Fse);
        Fmax = max(Fse)-Fmin;
    end
    
    Frel = (Fse-Fmin)/Fmax;
    
    if i == 1
        x1 = t(find(Frel > .9, 1));
    else
        x2 = t(Frel == max(Frel));
        y2 = max(Frel);
        x3 = max(t(Frel > .5*y2));
    end
    
    Fi(i,:) = interp1(t, Frel, ti);
    
end

X = [x2 y2 x3];

cost = sum((X - Y).^2);

figure(100)
plot(ti, Fi, '-', X(1), X(2), 'x', X(3), .5*X(2), 'x', ...
                  Y(1), Y(2), 'o', Y(3), .5*Y(2), 'o');  

drawnow
    
end
