clear all; close all; clc

[codefolder, name, ext] = fileparts(which('main_cfxc.m')); % folder contains the present m-file

remove_text = 0;
save_fig = 1;

figure(1)
color = get(gca,'colororder');

% can be used to generate either Fig 5 or Fig S2
figname = 'Fig5';
% figname = 'FigS2';

%% Human data
data = readmatrix('tetanus.csv');

Fmax_exp = max(data(:,2));
Frel_exp = data(:,2)/Fmax_exp;
texp = data(:,1)-.055;

subplot(242); plot(texp,Frel_exp,'-','linewidth',2,'color',[.8 .8 .8])
subplot(244); plot(texp,Frel_exp,'-','linewidth',2,'color',[.8 .8 .8])

subplot(246); plot(texp,Frel_exp,'--','linewidth',1,'color',[.5 .5 .5]); hold on
subplot(248); plot(texp,Frel_exp,'--','linewidth',1,'color',[.5 .5 .5]); hold on
% 
Frise(1:2,1) = min(texp(Frel_exp>.9));

data = readmatrix('twitch.csv');
texp_twitch = data(:,1);
Frel_exp_twitch = data(:,2)/Fmax_exp;
subplot(246); plot(data(:,1),Frel_exp_twitch,'-','linewidth',2,'color',[.8 .8 .8]); 
subplot(248); plot(data(:,1),Frel_exp_twitch,'-','linewidth',2,'color',[.8 .8 .8]); 

Fpeak(1:2,1) = max(Frel_exp_twitch);
tpeak(1:2,1) = data(data(:,2)==max(data(:,2)),1);

% make nice
figure(1)
for i = 1:8
    subplot(2,4,i); axis([0 .3 0 1]); box off; hold on
    
    xlabel('Time (s)')
    
    if mod(i,2),  ylabel('Activation')
    else, ylabel('Force')
    end

    if i == 6 || i == 8
        ylim([0 .5])
    end
end

%% models
parms.A = 1;
types = {'fast','slow'};

if strcmp(figname, 'Fig5') % compare Hill and Huxley
    models = {'Hill','Huxley'};
elseif strcmp(figname, 'FigS2') % compare Huxley to DM
    models = {'CB','Huxley'};
end

nmbs = [(1:2:13);
        (2:2:14)];

lt = {'--',':'};
lw = [2 1];

evs = 2:2:8;
ods = 1:2:9;

for m = 1:2
    for i = 1:2


        figure(1)
        load(['Fig5_',types{i},models{m},'_tetanus.mat'],'X','t');
        Fnet = X(:,2) - X(1,2);
        Fmax = max(Fnet);
        Frel = Fnet / Fmax;

        subplot(2,4,evs(i)); hold on; box off
        plot(t, Frel,lt{m},'color',color(1,:), 'linewidth',lw(m)); hold on
        
        subplot(2,4,ods(i)); hold on; box off
        plot(t, X(:,1)/parms.A,lt{m},'color',color(1,:), 'linewidth',lw(m)); hold on

        % stats
        Crise(m,i+1) = min(t(X(:,1)>.9*max(X(:,1))));
        Frise(m,i+1) = min(t(Frel>.9));
        Fint = interp1(t, Frel, texp);
        R2F(i,1,m) = 1 - sum((Fint-Frel_exp).^2,[],'omitnan') / sum((Frel_exp - mean(Frel_exp)).^2,[],'omitnan');

        figure(1)
        load(['Fig5_',types{i},models{m},'_twitch.mat'],'X','t');
        Fnet = X(:,2) - X(1,2);
        Frel = Fnet / Fmax;

        tpeak(m,i+1) = t(Frel(t<0.3)==max(Frel(t<0.3)));
        Fpeak(m,i+1) = max(Frel(t<0.3));
        Cpeak(m,i+1) = max(X(:,1));
        Fint = interp1(t, Frel, texp_twitch);
        R2F(i,2,m) = 1 - sum((Fint-Frel_exp_twitch).^2,[],'omitnan') / sum((Frel_exp_twitch - mean(Frel_exp_twitch)).^2,[],'omitnan');
        
        subplot(2,4,evs(i+2)); hold on; box off
        plot(t, Frel,lt{m},'color',color(1,:), 'linewidth',lw(m)); hold on
        
        subplot(2,4,ods(i+2)); hold on; box off
        plot(t, X(:,1)/parms.A, lt{m},'color',color(1,:), 'linewidth',lw(m)); hold on
%         ylim([0 .57])

    end
        model_Fpeak(m,:) = [Fpeak(m,2)/Fpeak(m,1) Fpeak(m,3)/Fpeak(m,1)];    
end


%% plot some dots
figure(1)
subplot(242); plot(Frise(1), .9, 'k.','markersize',10)
plot(Frise(1,2), .9, 'b.','markersize',10)
plot(Frise(1,1:2), [.8 .8], 'k-')
plot([Frise(1,1) Frise(1,1)], [.77 .83], 'k-')
plot([Frise(1,2) Frise(1,2)], [.77 .83], 'k-')

subplot(244); plot(Frise(1), .9, 'k.','markersize',10)
plot(Frise(1,3), .9, 'b.','markersize',10)
plot(Frise(1,[1,3]), [.8 .8], 'k-')

subplot(241); plot(Crise(1,2), .9, 'b.','markersize',10)
subplot(243); plot(Crise(1,3), .9, 'b.','markersize',10)

subplot(246); 
plot(tpeak(1,1), Fpeak(1,1), 'k.','markersize',10)
plot(tpeak(1,2), Fpeak(1,2), 'b.','markersize',10)
plot([tpeak(1,1) tpeak(1,1)], Fpeak(1,1:2),'k-')
plot(tpeak(1,1) + [-.007 .007], [Fpeak(1,2) Fpeak(1,2)],'k-')
plot([tpeak(1,1) tpeak(1,2)], [Fpeak(1,2) Fpeak(1,2)],'k:')

subplot(248); 
plot(tpeak(1,1), Fpeak(1,1), 'k.','markersize',10)
plot(tpeak(1,3), Fpeak(1,3), 'b.','markersize',10)
plot([tpeak(1,1) tpeak(1,1)], Fpeak(1,[1,3]),'k-')
plot(tpeak(1,1) + [-.007 .007], [Fpeak(1,3) Fpeak(1,3)],'k-')
plot([tpeak(1,1) tpeak(1,3)], [Fpeak(1,3) Fpeak(1,3)],'k:')

%% make nice
figure(1)
set(gcf,'position', [300 100 900 300])

for i = 1:8
    subplot(2,4,i); box off; hold on

    xlim([0 .3])
    set(gca,'xtick',0:.3:.3);
    
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
    saveas(gcf,figname,'jpg')
    disp(['Saved: ', [figname,'.jpg']])
end




