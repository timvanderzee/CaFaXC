clear all; close all; clc
saveFig = 1;

%% Code block: set Path.
%%% set the codefolder variable which is all the pathing that needs to be done.
% make sure that you're in the same folder as the main_cfxc.m file
[codefolder, name, ext] = fileparts(which('main_cfxc.m')); % folder contains the present m-file
% which we expect to be the current working directory.

savedir = fullfile(codefolder,'manuscript figures','model output');

% add sub-folders to path that contain data:
% 'kinematics_kinetics.mat' and 'muscle_excitations.mat'
addpath(fullfile(codefolder));
addpath(fullfile(codefolder,'model parameters'));
addpath(fullfile(codefolder,'cyclic force study','inputs'));

%% Update parameters
setts.mouse = 1; % false for human

if setts.mouse
    load('mouse_parms.mat')
    parms.exp.phi = 120; % joint angle (full extension = 0)
else
    parms.exp.phi = 90; % joint angle (full extension = 0)
    load('quad_parms.mat','parms','fl','fv');
end

% update
parms.ce.amin = 1e-3; % minimal excitation
parms = cfxc.calc_x0(parms); 

% needs to be stored here, because will be changed in next section
tau_original = parms.ce.tau;

%% Twitch and tetanus (Figs 4, 5, 7)
% parameters
parms.set.odeopt = odeset('maxstep',1e-3);
parms.set.no_tendon = 0;

parms.exp.vmtc = 0;
parms.exp.stim_type = 'u_func';
parms.exp.A = 1;
parms.exp.tstop = .3; 
parms.set.fixed_velocity = 0; % 

% settings
inp = [1 1 0 0];
sla = [0 1 0 1];

fignames = {'Tetanus - Calcium-based', 'Tetanus - Force-based','Twitch - Calcium-based', 'Twitch - Force-based'};
mtypes = {'Hill-type','crossbridge', 'CaFaXC'}; % if DM for crossbridge
mtypes{4} = 'Huxley'; % if you want to evaluate original model

titles = {'Ca activation', 'Force', 'CE length'};
ylabels = {'Activation', 'Force (N)', 'Length (m)'};

ss = 1:length(fignames);
ms = 1:length(mtypes);

for s = ss
    
    if ishandle(s), close(s); end 
    figure(s)
    set(gcf,'name',fignames{s});
    
    color = get(gca,'colororder');
    subplot(241); plot(fl.Lses, fl.Fses,'--','linewidth',1,'color',[.5 .5 .5]); hold on
    title('SEE force-length'); 
    ylim([0 parms.ce.Fmax])

    subplot(242); plot(fl.Lpes, fl.Fpes,'--','linewidth',1,'color',[.5 .5 .5]); hold on
    title('PE force-length'); 
    ylim([0 parms.ce.Fmax])

    subplot(243); plot(fl.Lces(:,2), fl.Fces(:,2),'--','linewidth',1,'color',[.5 .5 .5]); hold on
    title('CE force-length'); 
    ylim([0 parms.ce.Fmax])

    subplot(244); plot(fv.vHill, fv.FHill,'--','linewidth',1,'color',[.5 .5 .5]); hold on
    title('CE force-velocity'); 
    ylim([0 parms.ce.Fasymp])
    
    setts.step_input = inp(s);
    setts.slow_act = sla(s); % slow activation time constants

    fs = [linspace(.1,100,10) 200];
    Fpeak = ones(length(fs), 3);

    for c = 1 %:length(fs)
        parms.exp.freq = fs(c);

        if setts.step_input 
              parms.exp.u_func = @(t, parms) parms.exp.A .* (.5 + .5 * square(2*pi*parms.exp.freq*t)) .* (t < parms.exp.tstop);
        else, parms.exp.u_func = @(t, parms) parms.exp.A .* (.5 + .5 * square(2*pi*parms.exp.freq*t, .5 * parms.exp.freq)) .* (t < parms.exp.tstop); % assuming a pulse width of 5 ms
        end

        for m = 1:3

            % model type
            setts.M = m;

            % specify type
            parms.type = mtypes{m};
            
            if strcmp(parms.type,'Hill-type') % sometimes needed
                parms.set.odeopt = odeset('maxstep',1e-3);
            else
                parms.set.odeopt = odeset('maxstep',1e-3);
            end
            
%             parms.set.odeopt = [];
            
            if setts.slow_act && ~strcmp(parms.type,'CaFaXC')
               if setts.mouse
            %         parms.ce.tau = [.04 .012]; % fitted on force
                    parms.ce.tau(1) = .04;
               else
            %         parms.ce.tau = [.06 .045]; % fitted on force
                    parms.ce.tau(1) = .06;
               end
            else
                parms.ce.tau = tau_original;
            end

            if strcmp(parms.type,'crossbridge')
                X0 = parms.exp.x0;
            elseif  strcmp(parms.type,'CaFaXC')
                X0 = [parms.exp.x0(1) parms.exp.x0(1) parms.exp.x0(2:end)];
            elseif  strcmp(parms.type,'Hill-type')
                X0 = [parms.ce.amin parms.exp.l0]; 
            elseif strcmp(parms.type,'Huxley')
                X0 = [parms.exp.x0(1) zeros(1,parms.CB.nbins) parms.exp.x0(end)];
            end

            % simulate
            disp([fignames{s}, ' ', parms.type, ' model ']);
            tic
            [t,x] = ode113(@cfxc.sim_muscle, [0 parms.exp.tstop], X0, parms.set.odeopt, parms);
            d = toc;
            
            [y,X] = cfxc.get_model_output(t, x, parms);
       
            Fpeak(c,m) = max(y.Fse) - min(y.Fse);

            figure(s)
            subplot(241); plot(y.lse, y.Fse,'.','color',color(m,:)); 
            subplot(242); plot(y.lce, y.Fpe,'.','color',color(m,:)); 
            subplot(243); plot(y.lce, parms.func.fce(y.lce,parms)*parms.ce.Fmax,'.','color',color(m,:)); 

            for i = 1:3
                subplot(2,4,i+4); hold on
                plot(t,X(:,i),'color',color(m,:)); 
                title(titles{i});
                xlabel('Time (s)'); 
                ylabel(ylabels{i}); 
            end
            
            subplot(244)
            plot(y.vce/parms.ce.lceopt, y.Fce/parms.ce.Fmax, 'color',color(m,:));
            
            subplot(248);
            plot(t, y.vce, 'color',color(m,:)); 
            title('CE velocity'); hold on
            xlabel('Time (s)'); 
            ylabel('Velocity (m/s)');

            subplot(246);
            plot(t, y.Fpe,'--','color',color(m,:));
            plot(t, y.Fpe + y.Fce,'k:')
            ylim([0 parms.ce.Fmax])

            Fserel = (y.Fse - y.Fse(1)) / max((y.Fse - y.Fse(1)));
            
            
            if ~(strcmp(parms.type,'CaFaXC') && (s == 2 || s == 4)) % force-based force facilitation doesn't make sense
                
                if s < 3 % twitch
                    t90 = min(t(Fserel > .9));
                    disp(['90% rise time = ', num2str(t90*1000,3), ' ms'])
                
                else
                    tpeak = t(Fserel == max(Fserel));
                    disp(['Time to peak = ', num2str(tpeak*1000,3), ' ms'])
                end
                
            end
            
        disp(['Simulation time: ', num2str(round(d, 3)), ' s - Integration steps: ', num2str(length(t))])
            
        % saving
        save_as_Fig(saveFig, savedir, t, X, setts, parms);
        disp(' ')
        end
        

    end
    
    for i = 1:8
        subplot(2,4,i); box off;
    end
end


%%
function[] = save_as_Fig(saveFig, dir, tr, Xr, set, parms)

if saveFig
    
    
    Mnames = {'Hill', 'CB', '','Huxley'};
    
    
    if set.slow_act
        actname = 'slow';
    else
        actname = 'fast';
    end
    
    if set.mouse
        figname = 'Fig4';
    else
        figname = 'Fig5';
    end
    
    if set.M == 3
        if ~set.slow_act
            figname = 'Fig7';
            if set.mouse
                actname = 'mouse';
            else
                actname = 'human';
            end
        else
            return
        end
    end

    if set.step_input
        condname = 'tetanus';
    else
        condname = 'twitch';
    end

    % downsample
    t = linspace(0,max(tr),400);
    X = nan(length(t),size(Xr,2));
    for i = 1:size(Xr,2)
        X(:,i) = interp1(tr(:), Xr(:,i), t(:));
    end
    
    Fse = X(:,2);
    savename = [figname, '_', actname, Mnames{set.M}, '_',condname,'.mat'];
    fullsavename = fullfile(dir,figname,savename);
    save(fullsavename, 't','X','Fse','parms','set');
    
    disp(['Saved: ', fullsavename])
    
end
end
