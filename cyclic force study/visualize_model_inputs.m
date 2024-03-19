clear all; close all; clc
load('kinematics_kinetics.mat')
load('muscle_excitations.mat')
load('calcium_activations.mat')
load('crossbridge_activations.mat')

figure(1); 
color = [linspace(0,1,length(freqs))', zeros(length(freqs),2)];
for i = 1:length(freqs)
    
    subplot(221); plot(tx(i,:), U_opt_smooth(i,:),'color',color(i,:)); hold on
    subplot(222); plot(tx(i,:), C_opt_smooth(i,:),'color',color(i,:)); hold on
    subplot(223); plot(tx(i,:), R_opt(i,:),'color',color(i,:)); hold on
    subplot(224); plot(tx(i,:), X_opt(i,:),'color',color(i,:)); hold on
end

titles = {'Excitation','Calcium','Facilitation','Force'};
for i = 1:4
    subplot(2,2,i);
    box off
    title(titles{i});
    xlabel('Time (s)')
end

%% Time delay
for i = 1:length(freqs)
    Utau90 = find( U_opt_smooth(i,:) > (.9 * max(U_opt_smooth(i,:))),1);
    Ftau90 = find(  Fx(:) > (.9 * max(Fx(:))),1);
    tdelay(i) = tx(i,Ftau90) - tx(i,Utau90);
end

%% Gain vs. frequency
P = tdelay .* freqs * 360;
G = 1./max(U_opt_smooth,[],2);

figure(2)
subplot(211);
semilogx(freqs, G/G(1))
ylim([0 1.5]); box off

subplot(212);
semilogx(freqs, -P)
ylim([-150 0]); box off

titles = {'Gain','Phase'};
units = {'','deg'};

for i = 1:2
    subplot(2,1,i);
    box off
    title(titles{i});
    xlabel('Frequency (Hz)')
    ylabel([titles{i}, ' (',units{i},')']);
end