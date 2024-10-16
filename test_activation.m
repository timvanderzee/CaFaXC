close all; clc; clear all


parms.a = 1;

parms.f = [20 300];
parms.g = [20 300];
parms.u = 0.5;

parms.k = 50;


lw = {'-','--'};

ks = linspace(0,1,10);
color = [linspace(0,1,length(ks))' linspace(1,0,length(ks))', zeros(length(ks),1)];

parms.a = 1;
parms.u = .5;
tmax = .5;
% ks = 1;
parms.r = 1;

for i = 1:length(ks)
%     parms.r = ks(i);
parms.a = ks(i);

    parms.f(2) = 300;
    [t0,x0] = ode113(@dxfunc, [0 tmax], [0 0], [], parms);

    % set Q0 back to 0
    parms.f(2) = 0;
    [t1,x1] = ode113(@dxfunc, [0 tmax], [x0(end,1) 0], [], parms);

    % allow attachment
    parms.f(2) = 300;
    [t2,x2] = ode113(@dxfunc, [0 tmax], x1(end,:), [], parms);  
  
    % set a to max
    parms.a = 1;
    [t3,x3] = ode113(@dxfunc, [0 tmax], x2(end,:), [], parms);  

    t = [t0; t1+t0(end); t2+t0(end)+t1(end); t3+t0(end)+t1(end)+t2(end)];
    x = [x0; x1; x2; x3];

%     ttp(i) = min(t(x(:,2) > .9 * max(x(:,2))));

    figure(1)
    subplot(121)
    plot(t, x(:,1),'color',color(i,:)); hold on

    subplot(122)
    plot(t, x(:,2),'color',color(i,:)); hold on
end
% % 
% figure(2)
% plot(ks, ttp)

function[dx] = dxfunc(t, x, parms)


    N = x(1);
    Q = x(2);

    
    % Campbell 2018
%     Jon =  parms.f(1)/2 * (parms.r-N)   * (1 + parms.k * N);
%     Joff = parms.g(1)   * (N-Q/parms.a) * (1 + parms.k*(parms.r-N));
%     dQ = parms.f(2) * (parms.a*N-Q) - parms.g(2) * Q;
    
    Ntot = parms.a * parms.r;
    Jon =  parms.f(1)/2 * (Ntot-N)   * (1 + parms.k * N);
    Joff = parms.g(1)   * (N-Q) * (1 + parms.k*(Ntot-N));
    dQ = parms.f(2) * (N-Q) - parms.g(2) * Q;

    dN = Jon - Joff;
    dx = [dN; dQ];

end