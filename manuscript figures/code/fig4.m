clear all; close all; clc

[codefolder, name, ext] = fileparts(which('main_cfxc.m')); % folder contains the present m-file

remove_text = 0;
save_fig = 1;

figure(1)
color = get(gca,'colororder');

% can be used to generate either Fig 4 or Fig S1
figname = 'Fig4';
% figname = 'FigS1';

%% First-order
% showed first-order response as comparison
tlin = linspace(0,.2,1000);
tauC = .002;
tauF = .05;

Cfo = 1-exp(-tlin/tauC);
Ffo = 1-exp(-tlin/tauF);

ls = '--';

subplot(241); plot(tlin, Cfo,ls,'color',[.5 .5 .5], 'linewidth',1); 
subplot(242); plot(tlin, Ffo,ls,'color',[.5 .5 .5], 'linewidth',1);
subplot(243); plot(tlin, Cfo,ls,'color',[.5 .5 .5], 'linewidth',1); 
subplot(244); plot(tlin, Ffo,ls,'color',[.5 .5 .5], 'linewidth',1);
% subplot(246); plot(tlin, Ffo,'--','color',[.5 .5 .5], 'linewidth',1);
% subplot(248); plot(tlin, Ffo,'--','color',[.5 .5 .5], 'linewidth',1);

for i = 1:8
    subplot(2,4,i); 
    axis([0 .12 0 1.1]); hold on; box off
    xlabel('Time (s)')
    
    if mod(i,2),  ylabel('Activation')
    else, ylabel('Force')
    end

end

%% Data - tetanus
Ftetanus = readmatrix('H1996_Fig6_Force.csv');
Ctetanus = readmatrix('H1996_Fig6_CaTrop.csv');
Ftetanus(:,2) = (Ftetanus(:,2)-Ftetanus(1,2))/.27;
Ctetanus(:,2) = (Ctetanus(:,2)-Ctetanus(1,2))/.3;

Ftetanus(Ftetanus(:,1)>100,2) = nan;
Ctetanus(Ctetanus(:,1)>100,2) = nan;

% plot data
subplot(244); plot(Ftetanus(:,1)/1000,Ftetanus(:,2),'color',[.8 .8 .8], 'linewidth',2);
subplot(242); plot(Ftetanus(:,1)/1000,Ftetanus(:,2),'color',[.8 .8 .8], 'linewidth',2);
% subplot(243); plot(Ctetanus(:,1)/1000,Ctetanus(:,2),'color',[.8 .8 .8], 'linewidth',2);
% subplot(241); plot(Ctetanus(:,1)/1000,Ctetanus(:,2),'color',[.8 .8 .8], 'linewidth',2);

Cint = interp1(tlin, Cfo, Ctetanus(:,1)/1000);
R2C_for = 1 - sum((Cint-Ctetanus(:,2)).^2) / sum((Ctetanus(:,2) - mean(Ctetanus(:,2)).^2));

tF100(1:2,1) = min(tlin(Ffo>.9));
tC100(1:2,1) = min(tlin(Cfo>.9));

%% Data - twitch
Ftwitch = readmatrix('H1996_Fig4_twitch_Force.csv');
F2 = readmatrix('H1996_Fig4_tetanus_Force.csv');

Ftwitch(:,2) = (Ftwitch(:,2)-Ftwitch(1,2))/.75;
F2(:,2) = (F2(:,2)-F2(1,2))/.75;

subplot(248); plot(Ftwitch(:,1)/1000,Ftwitch(:,2),'color',[.8 .8 .8], 'linewidth',2); ylim([0 .5])
subplot(246); plot(Ftwitch(:,1)/1000,Ftwitch(:,2),'color',[.8 .8 .8], 'linewidth',2); ylim([0 .5])
subplot(247); plot(Ctetanus(1:30,1)/1000, Ctetanus(1:30,2),'color',[.8 .8 .8], 'linewidth',2);
subplot(245); plot(Ctetanus(1:30,1)/1000, Ctetanus(1:30,2),'color',[.8 .8 .8], 'linewidth',2);

Fpeak(1:2,1) = max(Ftwitch(:,2));
[Cpeak(1),i] = max(Ctetanus(1:30,2));

tFpeak(1:2,1) = Ftwitch((Ftwitch(:,2) == max(Ftwitch(:,2))),1)/1000;
tCpeak(1:2,1) = Ctetanus((Ctetanus(1:30,2) == max(Ctetanus(1:30,2))),1)/1000;

tauC = .012;
tlin = tlin(:);
Cdecay = exp(-(tlin-Ctetanus(i)/1000)/tauC);

subplot(247); plot(tlin, Cdecay,'--','color',[.5 .5 .5], 'linewidth',1);
subplot(245); plot(tlin, Cdecay,'--','color',[.5 .5 .5], 'linewidth',1);

% combine decay and original
tcomb = [Ctetanus(1:30,1)/1000; tlin(tlin > (Ctetanus(30,1)/1000))];
Ccomb = [Ctetanus(1:30,2); Cdecay(tlin>(Ctetanus(30,1)/1000))];

%% Model simulations - Tetanus
figure(1)

if strcmp(figname, 'Fig4') % compare Hill and Huxley
    models = {'Hill','Huxley'};
elseif strcmp(figname, 'FigS1') % compare Huxley to DM
    models = {'CB','Huxley'};
end

lt = {'--',':'};
% lt = {'-','-'};
lw = [2 1];

for m = 1:2
    model = models{m};

    load(['Fig4_slow',model,'_tetanus.mat'],'Fse','t','X')
    
    Fmax = max(Fse);
    Fmin = Fse(1);
    Frel = (Fse-Fmin) / (Fmax-Fmin);

    subplot(244); plot(t, Frel,lt{m},'color',color(1,:), 'linewidth',lw(m));
    subplot(243); plot(t, X(:,1),lt{m},'color',color(1,:), 'linewidth',lw(m));

    % stats
    Fint = interp1(t, Frel, Ftetanus(:,1)/1000);
    Cint = interp1(t, X(:,1), Ctetanus(:,1)/1000);
    % subplot(244); plot(F(:,1)/1000, Fint, '.')
    R2F(1,1,m) = 1 - sum((Fint-Ftetanus(:,2)).^2,[],'omitnan') / sum((Ftetanus(:,2) - mean(Ftetanus(:,2),'omitnan').^2),[],'omitnan');
    R2C(1,1,m) = 1 - sum((Cint-Ctetanus(:,2)).^2,[],'omitnan') / sum((Ctetanus(:,2) - mean(Ctetanus(:,2),'omitnan').^2),[],'omitnan');
    tF100(m,2) = min(t(Frel > .9));
    tC100(m,2) = min(t(X(:,1) > .9));

    figure(1)
    load(['Fig4_fast',model,'_tetanus.mat'],'Fse','t','X')
    Fmax = max(Fse);
    Frel = (Fse-Fmin) / (Fmax-Fmin);
    subplot(242); plot(t, Frel,lt{m},'color',color(1,:), 'linewidth',lw(m));
    subplot(241); plot(t, X(:,1),lt{m},'color',color(1,:), 'linewidth',lw(m));

    Fint = interp1(t, Frel, Ftetanus(:,1)/1000);
    Cint = interp1(t, X(:,1), Ctetanus(:,1)/1000);
    % subplot(244); plot(F(:,1)/1000, Fint, '.')
    R2F(1,2,m) = 1 - sum((Fint-Ftetanus(:,2)).^2,[],'omitnan') / sum((Ftetanus(:,2) - mean(Ftetanus(:,2),'omitnan').^2),[],'omitnan');
    R2C(1,2,m) = 1 - sum((Cint-Ctetanus(:,2)).^2,[],'omitnan') / sum((Ctetanus(:,2) - mean(Ctetanus(:,2),'omitnan').^2),[],'omitnan');

    tF100(m,3) = min(t(Frel > .9));
    tC100(m,3) = min(t(X(:,1) > .9));

    % Quantify lead
    model_Flead(m,:) = [tF100(m,1)-tF100(m,2) tF100(m,1)-tF100(m,3)];
    model_Clead(m,:) = [tC100(m,1)-tC100(m,2) tC100(m,1)-tC100(m,3)];


%% model simulations - Twitch
    model = models{m};
    load(['Fig4_slow',model,'_twitch.mat'],'Fse','t','X')
    Frel = (Fse-Fmin) / (Fmax-Fmin);
    subplot(248); plot(t, Frel,lt{m},'color',color(1,:), 'linewidth',lw(m)); %ylim([0 .3])
    subplot(247); plot(t, X(:,1),lt{m},'color',color(1,:), 'linewidth',lw(m));

    Fint = interp1(t, Frel, Ftwitch(:,1)/1000);
    Cint = interp1(t, X(:,1), tcomb);
    % subplot(244); plot(F(:,1)/1000, Fint, '.')
    R2F(2,1,m) = 1 - sum((Fint-Ftwitch(:,2)).^2,[],'omitnan') / sum((Ftwitch(:,2) - mean(Ftwitch(:,2))).^2,[],'omitnan');
    R2C(2,1,m) = 1 - sum((Cint-Ccomb).^2,[],'omitnan') / sum((Ccomb - mean(Ccomb)).^2,[],'omitnan');

    tFpeak(m,2) = t(Frel(t<0.1) == max(Frel(t<0.1)));
    tCpeak(m,2) = t(X(:,1) == max(X(:,1)));
    Fpeak(m,2) = max(Frel(t<0.1));
    Cpeak(m,2) = max(X(:,1));

    load(['Fig4_fast',model,'_twitch.mat'],'Fse','t','X')
    Frel = (Fse-Fmin) / (Fmax-Fmin);
    subplot(246); plot(t, Frel,lt{m},'color',color(1,:), 'linewidth',lw(m)); %ylim([0 .3])
    subplot(245); plot(t, X(:,1),lt{m},'color',color(1,:), 'linewidth',lw(m));

    Fint = interp1(t, Frel, Ftwitch(:,1)/1000);
    Cint = interp1(t, X(:,1), tcomb);
    % subplot(244); plot(F(:,1)/1000, Fint, '.')
    R2F(2,2,m) = 1 - sum((Fint-Ftwitch(:,2)).^2,[],'omitnan') / sum((Ftwitch(:,2) - mean(Ftwitch(:,2))).^2,[],'omitnan');
    R2C(2,2,m) = 1 - sum((Cint-Ccomb).^2,[],'omitnan') / sum((Ccomb - mean(Ccomb)).^2,[],'omitnan');

    tFpeak(m,3) = t(Frel(t<0.1) == max(Frel(t<0.1)));
    Fpeak(m,3) = max(Frel(t<0.1));
    tCpeak(m,3) = t(X(:,1) == max(X(:,1)));
    Cpeak(m,3) = max(X(:,1));

    % Quantify peak force
    model_Fpeak(m,:) = [Fpeak(m,2)/Fpeak(m,1) Fpeak(m,3)/Fpeak(m,1)];
    mode_Cpeak(m,:) = [Cpeak(m,2)/Cpeak(m,1) Cpeak(m,3)/Cpeak(m,1)];
end

%% dots and lines
figure(1)
subplot(244); plot(tF100(1,1), .9, 'k.','markersize',10)
plot(tF100(1,2), .9, 'b.','markersize',10)
plot(tF100(1,1:2), [.8 .8], 'k-')

subplot(242); plot(tF100(1,1), .9, 'k.','markersize',10)
plot(tF100(1,3), .9, 'b.','markersize',10)
plot(tF100(1,[1,3]), [.8 .8], 'k-')

subplot(243); plot(tC100(1,1), .9, 'k.','markersize',10)
plot(tC100(1,2), .9, 'b.','markersize',10)
plot(tC100(1,1:2), [.8 .8], 'k-')

subplot(241); plot(tC100(1), .9, 'k.','markersize',10)
plot(tC100(1,3), .9, 'b.','markersize',10)
plot(tC100(1,[1,3]), [.8 .8], 'k-')

subplot(248); 
plot(tFpeak(1,1), Fpeak(1,1), 'k.','markersize',10)
plot(tFpeak(1,2), Fpeak(1,2), 'b.','markersize',10)
plot([tFpeak(1,1) tFpeak(1,1)], Fpeak(1,1:2),'k-')
plot([tFpeak(1,1) tFpeak(1,2)], [Fpeak(1,2) Fpeak(1,2)],'k--')

subplot(246); 
plot(tFpeak(1,1), Fpeak(1,1), 'k.','markersize',10)
plot(tFpeak(1,3), Fpeak(1,3), 'b.','markersize',10)
plot([tFpeak(1,1) tFpeak(1,1)], Fpeak(1,[1,3]),'k-')
plot([tFpeak(1,1) tFpeak(1,2)], [Fpeak(1,3) Fpeak(1,3)],'k--')

subplot(247); 
plot(tCpeak(1,1), Cpeak(1,1), 'k.','markersize',10)
plot(tCpeak(1,2), Cpeak(1,2), 'b.','markersize',10)
plot([tCpeak(1,1) tCpeak(1,1)]+.01, Cpeak(1,1:2),'k-')
plot([tCpeak(1,1)+.01 tCpeak(1,2)], [Cpeak(1,2) Cpeak(1,2)],'k--')
plot([tCpeak(1,1)+.01 tCpeak(1,2)], [Cpeak(1,1) Cpeak(1,1)],'k--')

subplot(245); 
plot(tCpeak(1,1), Cpeak(1,1), 'k.','markersize',10)
plot(tCpeak(1,3), Cpeak(1,3), 'b.','markersize',10)
plot([tCpeak(1,1) tCpeak(1,1)], Cpeak(1,[1,3]),'k-')
plot([tCpeak(1,1) tCpeak(1,3)], [Cpeak(1,3) Cpeak(1,3)],'k--')

%% make nice
figure(1)
set(gcf,'position', [300 100 900 300])

for i = 1:8
    subplot(2,4,i);
    set(gca,'xtick',[0 .1]);
    set(gca,'ytick',0:.5:1);
    
    if remove_text
        xlabel(''); ylabel('')
        set(gca,'xticklabels',[])
        set(gca,'yticklabels',[])
    end
end

%% save
if save_fig
    cd(codefolder)
    saveas(gcf,figname,'jpg')
    disp(['Saved: ', [figname,'.jpg']])
end

