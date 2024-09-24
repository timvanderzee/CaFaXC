clear all; close all; clc

[codefolder, name, ext] = fileparts(which('main_cfxc.m')); % folder contains the present m-file

remove_text = 0;
save_fig = 1;

figure(1)
color = get(gca,'colororder');

%% Human data
dt = .01;

data = readmatrix('tetanus.csv');

Fmax_exp = max(data(:,2));
Frel_exp = data(:,2)/Fmax_exp;
texp1 = data(:,1)-.055 + dt;

subplot(244); 
plot(texp1,Frel_exp,'-','linewidth',2,'color',[.8 .8 .8]); hold on

subplot(248); 
plot(texp1,Frel_exp,'-','linewidth',1,'color',[.5 .5 .5]); hold on

Frise = min(texp1(Frel_exp>.9));

data = readmatrix('twitch.csv');
texp = data(:,1) + dt;

subplot(248); plot(texp,data(:,2)/Fmax_exp,'-','linewidth',2,'color',[.8 .8 .8]); ylim([0 .5])

tpeak(1,1) = texp(data(:,2)==max(data(:,2)));
Fpeak(1,1) = max(data(:,2)/Fmax_exp);

%% Mouse tetanus data
tlin = linspace(0,.2,1000);
tauC = .002;
tauF = .05;

Cfo = 1-exp(-tlin/tauC);
Ffo = 1-exp(-tlin/tauF);
Frise(1,2) = min(tlin(Ffo>.9));

F = readmatrix('H1996_Fig6_Force.csv');
C = readmatrix('H1996_Fig6_CaTrop.csv');
F(:,2) = (F(:,2)-F(1,2))/.27;
C(:,2) = (C(:,2)-C(1,2))/.3;

F(F(:,1)>100,2) = nan;
C(C(:,1)>100,2) = nan;

subplot(241); 
plot(tlin,Cfo,'--','color',[.5 .5 .5], 'linewidth',1);

subplot(242); 
plot(F(:,1)/1000,F(:,2),'color',[.8 .8 .8], 'linewidth',2); hold on
plot(tlin, Ffo,'--','color',[.5 .5 .5], 'linewidth',1); hold on

%% Mouse twitch data
Mtwitch = readmatrix('H1996_Fig4_twitch_Force.csv');
Mtwitch(:,2) = (Mtwitch(:,2)-Mtwitch(1,2))/.75;

tpeak(1,2) = Mtwitch(Mtwitch(:,2)==max(Mtwitch(:,2)),1)/1000;
Fpeak(1,2) = max(Mtwitch(:,2));

tlin = linspace(0,.1,1000);
tlin = tlin(:);

[Cpeak(1),i] = max(C(1:30,2));
tauC = .012;
Cdecay = exp(-(tlin-C(i)/1000)/tauC);

subplot(245); 
plot(tlin, Cdecay,'color',[.8 .8 .8], 'linewidth',2);

subplot(246); 
plot(Mtwitch(:,1)/1000,Mtwitch(:,2),'color',[.8 .8 .8], 'linewidth',2); hold on; ylim([0 .5])

%% Human simulations
parms.A = .57;
parms.A = 1;
load('Fig7_human_tetanus.mat','X','t')
Fmax = max(X(:,2))-min(X(:,2)); 
Cmax = max(X(:,1));

Frel = (X(:,2)-X(1,2)) / Fmax;
Crel = X(:,1) / Cmax;

% stats
Rrise(1,1) = min(t(X(:,4)>.9));
Frise(2,1) = min(t(Frel>.9));
Fint = interp1(t, Frel, texp1);
R2F(1,1) = 1 - sum((Fint-Frel_exp).^2,[],'omitnan') / sum((Frel_exp - mean(Frel_exp)).^2,[],'omitnan');

subplot(244);  
plot(t, Frel,'--','color',color(1,:),'linewidth',2); hold on;
plot(t,X(:,4)/parms.A,':','color',color(1,:),'linewidth',1)

subplot(243); 
plot(t, Crel,'--','color',color(1,:),'linewidth',2); hold on; 

load('Fig7_human_twitch.mat','X','t')
Frel_exp = data(:,2)/Fmax_exp;
Frel = (X(:,2)-X(1,2)) / Fmax;
Crel = X(:,1) / Cmax;

subplot(244);
plot(t, Frel,'--','color',color(1,:),'linewidth',2);  hold on;

subplot(248);
plot(t, Frel,'--','color',color(1,:),'linewidth',2);  hold on;
plot(t,X(:,4)/parms.A,':','color',color(1,:),'linewidth',1)

subplot(247); 
plot(t, Crel,'--','color',color(1,:),'linewidth',2); ylim([0 1]); hold on; 

% stats
Fint = interp1(t, Frel, data(:,1));
R2F(1,2) = 1 - sum((Fint-Frel_exp).^2,[],'omitnan') / sum((Frel_exp - mean(Frel_exp)).^2,[],'omitnan');
tpeak(2,1) = t(Frel==max(Frel));
Fpeak(2,1) = max(Frel);

%% Mouse simulations
parms.A = 1;
load('Fig7_mouse_tetanus.mat','X','t')
Fmax = max(X(:,2)) - min(X(:,2)); 
Cmax = max(X(:,1));
Frel = (X(:,2)-X(1,2)) / Fmax;
Crel = X(:,1) / Cmax;

subplot(241); 
plot(t, Crel,'--','color',color(1,:),'linewidth',2); hold on; 

subplot(242); 
plot(t, Frel,'--','color',color(1,:),'linewidth',2); hold on; 
plot(t,X(:,4)/parms.A,':','color',color(1,:),'linewidth',1)

% stats
Fint = interp1(t, Frel, F(:,1)/1000);
R2F(2,1) = 1 - sum((Fint-F(:,2)).^2,[],'omitnan') / sum((F(:,2) - mean(F(:,2),'omitnan')).^2,[],'omitnan');
Rrise(1,2) = min(t(X(:,4)>.9));
Frise(2,2) = min(t(Frel>.9));

Cint = interp1(t, Crel, C(:,1)/1000);
Cint((C(:,1)/1000) > .1) = nan;
R2C(1,1) = 1 - sum((Cint-C(:,2)).^2,[],'omitnan') / sum((C(:,2) - mean(C(:,2))).^2,[],'omitnan');

load('Fig7_mouse_twitch.mat','X','t')
Frel = (X(:,2)-X(1,2)) / Fmax;
Crel = X(:,1) / Cmax;

subplot(246); 
plot(t, Frel,'--','color',color(1,:),'linewidth',2); hold on; 
plot(t,X(:,4)/parms.A,':','color',color(1,:),'linewidth',1)

% combine decay and original
tcomb = [C(1:30,1)/1000; tlin(tlin > (C(30,1)/1000))];
Ccomb = [C(1:30,2); Cdecay(tlin>(C(30,1)/1000))];

subplot(245); 
plot(tcomb, Ccomb, '--', 'color', [.5 .5 .5]);
plot(t, Crel,'--','color',color(1,:),'linewidth',2); ylim([0 1]); hold on; 

% stats
Fint = interp1(t, Frel, Mtwitch(:,1)/1000);
R2F(2,2) = 1 - sum((Fint-Mtwitch(:,2)).^2,[],'omitnan') / sum((Mtwitch(:,2) - mean(Mtwitch(:,2))).^2,[],'omitnan');
tpeak(2,2) = t(Frel==max(Frel));
Fpeak(2,2) = max(Frel);

Cint = interp1(t, Crel, tcomb);
Cint((C(:,1)/1000) > .1) = nan;
R2C(1,2) = 1 - sum((Cint-Ccomb).^2,[],'omitnan') / sum((Ccomb - mean(Ccomb)).^2,[],'omitnan');

%% Make nice
figure(1)
for i = 1:8
    subplot(2,4,i); box off; hold on
    
    if mod(i,2) == 0
        xlim([0 .3])
    else
        xlim([0 .12])
    end
    
end



%% plot some dots
figure(1)
subplot(242); plot(Frise(1,2), .9, 'k.','markersize',10)
plot(Frise(2,2), .9, 'b.','markersize',10)
plot(Frise(1:2,2), [.8 .8], 'k-')
plot(Rrise(2), .9, 'g.')

subplot(244); plot(Frise(1,1), .9, 'k.','markersize',10)
plot(Frise(2,1), .9, 'b.','markersize',10)
plot(Frise(1:2,1), [.8 .8], 'k-')
plot(Rrise(1), .9, 'g.')
subplot(248); 
plot(tpeak(1,1), Fpeak(1,1), 'k.','markersize',10)
plot(tpeak(2,1), Fpeak(2,1), 'b.','markersize',10)

subplot(246); 
plot(tpeak(1,2), Fpeak(1,2), 'k.','markersize',10)
plot(tpeak(2,2), Fpeak(2,2), 'b.','markersize',10)

%% Make nice
figure(1)
set(gcf,'position', [300 100 900 300])

for i = 1:8
    subplot(2,4,i); box off; hold on

    % mouse vs. human
    if sum(i == [3 4 7 8])>0
        xlim([0 .3])
        set(gca,'xtick',0:.3:.3);
    else
        xlim([0 .12])
        set(gca,'xtick',0:.1:.3);
    end
    
    % activation vs. force
    if mod(i,2) == 0
        ylabel('Force')
    else
        ylabel('Activation')
    end
    
    if i < 5, set(gca,'xticklabels',[])
    else
        xlabel('Time (s)')
    end
    
    if i == 6 || i == 8
        ylim([0 .5])
        set(gca,'ytick',0:.5:1);
    else     
        ylim([0 1.05])
        set(gca,'ytick',0:1:1);
    end
    
    if remove_text
        xlabel(''); ylabel('')
        set(gca,'xticklabels',[])
        set(gca,'yticklabels',[])
    end
end

%% save
if save_fig
    cd(codefolder)
    saveas(gcf,'Fig7','jpg')
    disp('Saved: Fig7.jpg')
end