clear all; close all; clc

main = 'C:\Users\timvd\OneDrive\4. Calcium - TBD\new paper';
% main = 'C:\Users\tim.vanderzee\OneDrive\4. Calcium - TBD\new paper';

cd([main,'\PNG\data'])
F = readmatrix('BH1996_Force.csv');
C = readmatrix('BH1996_CaTrop.csv');

%
F(:,2) = (F(:,2)-F(1,2))/.27;
C(:,2) = (C(:,2)-C(1,2))/.3;

F(F(:,1)>100,2) = nan;
C(C(:,1)>100,2) = nan;

%% first-order
% showed first-order response as comparison
tlin = linspace(0,.2,1000);
tauC = .002;
tauF = .05;

Cfo = 1-exp(-tlin/tauC);
Ffo = 1-exp(-tlin/tauF);

subplot(221); plot(tlin, Cfo,'b--','linewidth',2); hold on; ylabel('Activation')
subplot(222); plot(tlin, Cfo,'b--','linewidth',2); hold on

id = find(Cfo>.9,1);
tc100 = tlin(id);
plot(tc100, .9,'k.','markersize',10)

id = find(Ffo>.9,1);
t100(1) = tlin(id);

plot([tc100 t100(1)], [.9 .9],'k-')

for i = 1:3
    subplot(2,2,i+1); plot(tlin, Ffo,'r:','linewidth',1); hold on
    plot(t100(1), .9,'k.','markersize',10)
end

%% data
subplot(221); plot(C(:,1)/1000,C(:,2),'color',[.6 .6 1], 'linewidth',2);

for i = 2:4
    subplot(2,2,i); plot(F(:,1)/1000,F(:,2),'color',[1 .6 .6],'linewidth',2);
end

titles = {'Activation dynamics', 'Force development dynamics', 'Hill-type model', 'Cross-bridge model'};

for i = 1:4
    subplot(2,2,i); axis([0 .12 0 1.1]); hold on; box off
    title(titles{i}); xlabel('Time')
    set(gca,'xtick',[])
    set(gca,'ytick',[0 1])
end

subplot(222); plot([.05 .07], [.3 .3], 'color',[.5 .5 .5], 'linewidth',5)



%% model simulations
cd([main,'\van der Zee & Kuo 2023\figures\MAT'])
load('Fig3_Hilltype.mat')
Fserel = (Fse - Fse(1)) / max((Fse - Fse(1)));

subplot(223); plot(t, Fserel,'r--','linewidth',2);

id = find(Fserel>.9,1);
t100(2) = t(id);

plot(t100(2), Fserel(id),'k.','markersize',10)

load('Fig3_crossbridge.mat')
Fserel = (Fse - Fse(1)) / max((Fse - Fse(1)));
subplot(224); plot(t, Fserel,'r--','linewidth',2);

id = find(Fserel>.9,1);
t100(3) = t(id);

plot(t100(3), Fserel(id),'k.','markersize',10)

%% Quantify lead
Hill_lead = t100(1)-t100(2)
CB_lead = t100(1)-t100(3)

subplot(223); plot([t100(1) t100(2)], [.9 .9],'k-'); ylabel('Force')
subplot(224); plot([t100(1) t100(3)], [.9 .9],'k-')







