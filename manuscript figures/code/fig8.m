clear all; close all; clc

save_fig = 1;

[codefolder, name, ext] = fileparts(which('main_cfxc.m')); % folder contains the present m-file

cd(codefolder)
addpath(genpath(codefolder))
load('quad_parms.mat')

figure(1)
color = get(gca,'colororder');

%% Twitch and tetanus
% same as Fig 7
tetanus = readmatrix('tetanus.csv');
twitch = readmatrix('twitch.csv');

Fmax_exp = max(tetanus(:,2));
Frel_exp = tetanus(:,2)/Fmax_exp;
texp = tetanus(:,1)-.055;

figure(1); 
subplot(322); 
plot(texp,Frel_exp,'-','linewidth',2,'color',[.5 .5 .5]); hold on
plot(twitch(:,1),twitch(:,2)/Fmax_exp,'-','linewidth',2,'color',[.8 .8 .8]); 
axis([0 .3 0 1]); box off; xlabel('Time (s)')
title('Twitch & tetanus')

load('Fig7_human_tetanus.mat','X','t')
Fmax = max(X(:,2))-min(X(:,2)); 
Frel = (X(:,2)-X(1,2)) / Fmax;

R2(1) = calc_R2(texp, Frel_exp, t, Frel);

subplot(322); hold on; 
plot(t, Frel,'-','color',color(1,:),'linewidth',2)

load('Fig7_human_twitch.mat','X','t')
Frel = (X(:,2)-X(1,2)) / Fmax;
plot(t, Frel,'-','color',color(1,:),'linewidth',2)

R2(2) = calc_R2(twitch(:,1), twitch(:,2)/Fmax_exp, t, Frel);



%% Force-velocity
figure(1)
TV = readmatrix('force_velocity.csv'); % torque-angular velocity
F = TV(:,2) / parms.r;
Fmax = F(4); % max isometric
v = TV(:,1)/180*pi * parms.r;

subplot(321); plot(v, F/Fmax,'o','color',[.5 .5 .5], 'markerfacecolor',[.5 .5 .5],'markersize',5); box off; hold on
xlabel('Velocity (m/s)');ylabel('Force'); title('Force-velocity')

load('Fig8_force-velocity.mat');
plot(-vall, Fss/parms.Fmax,'-','color',color(1,:),'linewidth',2)
xlim([-.3 .6])

R2(3) = calc_R2(v, F/Fmax, -vall, Fss/parms.Fmax); 

%% Force-length
FL = readmatrix('force_length.csv');

% re-calibrate (0 was actually 5)
L = (FL(:,1) / 12 * 7 + 5)/100;
Frel = FL(:,2)/max(FL(:,2));

subplot(323);
plot(L, Frel,'o','color',[.5 .5 .5], 'markerfacecolor',[.5 .5 .5],'markersize',5); hold on

load('Fig8_force-length.mat')
plot(Lmin, Fce_max/parms.Fmax,'-','color',color(1,:),'linewidth',2); box off
ylim([0 1]); xlabel('Length (m)'); ylabel('Force'); title('Force-length')

R2(4) = calc_R2(L, Frel, Lmin, Fce_max/parms.Fmax); 

%% Force-firing rate
FF = readmatrix('force_frequency.csv');

subplot(325); plot(FF(:,1), FF(:,2),'o','color',[.5 .5 .5], 'markerfacecolor',[.5 .5 .5],'markersize',5); box off; hold on
xlabel('Firing rate (Hz)');ylabel('Force')

load('Fig8_force-frequency.mat')
plot(freqs, Fpeak_freqs/(max(Fpeak_freqs)),'-','color',color(1,:),'linewidth',2)
axis([0 100.5 0 1])
title('Force - firing rate')

R2(5) = calc_R2(FF(:,1), FF(:,2), freqs, Fpeak_freqs/(max(Fpeak_freqs))); 

%% Low-pass filter
figure(1)
% data
efreqs = .5:.5:2.5;
load('lowpass_filter.mat');

Tmax70 = [ 1.2060    2.1846    1.7538    0.9421    1.8541      1.3585    1.4346 0.8822    2.1236] * 145;
peakemg = squeeze(mean(peakemg,1));
gain_exp_raw = (peaktrq./(Tmax70*2)) ./ peakemg;

% subtract individal offset
gain_exp = gain_exp_raw - mean(gain_exp_raw) + mean(gain_exp_raw(:));
phase_exp = -time_delay' * 360 .* efreqs;

subplot(324);   
semilogx(efreqs, mean(gain_exp,2),'o','color',[.5 .5 .5], 'markerfacecolor',[.5 .5 .5],'markersize',5); ylim([0 2.7]); hold on

subplot(326);
semilogx(efreqs, mean(phase_exp,1),'o','color',[.5 .5 .5], 'markerfacecolor',[.5 .5 .5],'markersize',5); ylim([-150 20]); hold on 

% model
load('Fig8_lowpass-filter','Fgain','Fphase','freqs')
 
% scale
scale = mean(gain_exp(:)) / mean(Fgain);
subplot(324); semilogx(freqs, Fgain*scale,'-','color',color(1,:),'linewidth',2); box off; 
 title('Low-pass filter gain'); xlabel('Contraction frequency (Hz)')

subplot(326); semilogx(freqs, Fphase,'-','color',color(1,:),'linewidth',2); box off
axis([.4 2.6 -150 0]); title('Low-pass filter phase'); xlabel('Contraction frequency (Hz)')

R2(6) = calc_R2(efreqs, mean(gain_exp,2), freqs, Fgain*scale); 
R2(7) = calc_R2(efreqs, mean(phase_exp,1), freqs, Fphase); 


%% make nice
for i = 1:5
    subplot(3,2,i)
    ylabel('Force');
end

subplot(326); ylabel('Phase (deg)')

%% figure size
set(gcf,'position', [0 0 400 450])

if save_fig
cd(codefolder)
saveas(gcf,'Fig8','jpg')
end

%% functions
function[R2] = calc_R2(x_exp, y_exp, x_sim, y_sim)
    
    % interp1 to experimental time points
    y_int = interp1(x_sim(:), y_sim(:), x_exp(:));
    
    % determine R2
    R2 = 1 - sum((y_exp(:)-y_int(:)).^2,[],'omitnan') / sum((y_exp(:) - mean(y_exp(:))).^2,[],'omitnan');
end


