classdef cfxc < handle
  properties (Access = public)
  end
  methods
    function obj=cfxc() %walk analysis
    end
  end

  methods(Static)            % ------STATIC------ %
    

    % moved from 'find_model_inputs'
    function[Uc, Xc] = solve_optimal_control(prob, parms, savedir)
    % given prob: -> freq, N, 
    % compute optimal Uc that optimally reproduces the desired frequency.

      for f = length(prob.freqs):-1:1
      % sine characteristics
      freq = prob.freqs(f);
      t = linspace(0,1/freq, prob.N);
      dt = mean(diff(t));
      
      opti = casadi.Opti(); 
      U = opti.variable(1,prob.N);
      X = opti.variable(prob.M,prob.N);
      
      disp('Dynamic contraints ...')
      
      for i = 1:prob.N-1
     
          Xk = X(:,i+0);
          Xk_plus = X(:,i+1);
      
          dXnow = cfxc.dXfunc(U(i+0), Xk, parms(i,f));
          dXnex = cfxc.dXfunc(U(i+1), Xk_plus, parms(i,f));
      
          Xhalf = 0.5*(Xk + Xk_plus) + dt/8 * (dXnow - dXnex);
          Uhalf = 0.5*(U(i) + U(i+1));
      
          dXhalf = cfxc.dXfunc(Uhalf, Xhalf, parms(i,f));
      
          xdot = (dXnow + 4*dXhalf + dXnex) / 6;
      
          opti.subject_to(cfxc.eulerIntegrator(Xk, Xk_plus, xdot, dt) == 0);
      
      end
      
       disp('Cost function ...')
      
      % cost function
      J = 0;
      
      if strcmp(parms(i,f).type,'crossbridge')
          Y = X(2,:)/parms(i,f).CB.Xmax(2); % relative force
      else
          Y = X;
      end
      
      for i = 1:prob.N
         J = J + (Y(i) - prob.Target(i,f)).^2;
      end
      
          disp('Initial guess ...')
      % initial conditions
          if prob.prev == 1
              if f == length(prob.freqs)
                  load(fullfile(savedir,'prev.mat'),'R_opt','X_opt');
                  load(fullfile(savedir,'crossbridge_activations'),'R_opt');
                  R_opt = R_opt(end,:);
              else
                  R_opt = U_opt;
              end
              
      
                for i = 1:prob.N
                  opti.set_initial(U(i), R_opt(i));
               end
              
              for j = 1:prob.M
                  for i = 1:prob.N
                      opti.set_initial(X(j,i), X_opt(j,i));
                  end
              end
          else
                for i = 1:prob.N
                  opti.set_initial(U(i), parms(i,f).Umin);
               end
              
              for j = 1:prob.M
                  for i = 1:prob.N
                      opti.set_initial(X(j,i), parms(i,f).Xmin(j));
                  end
              end
          end
          
          
          disp('Static constraints ...')
          
          for j = 1:prob.M
            opti.subject_to(X(j,1) == parms(i,f).Xmin(j));
          end
      
          % input constraints
          opti.subject_to(U(1) == parms(i,f).Umin);
          opti.subject_to(U >= parms(i,f).Umin);
          opti.subject_to(U <= 1);
      
          disp('Set cost function ...')
          % Set cost function.
          opti.minimize(J);
      
          p_opts = struct('expand',true);
          s_opts = struct('max_iter', 10000);
          opti.solver('ipopt',p_opts,s_opts);
      
          disp('Figure ...')
          figure(100); opti.callback(@(i) plot(opti.debug.value(Y)))
      
          disp('Solving ...')
          
          % Solve the NLP.
          try sol = opti.solve();
          catch

          end
      
          X_opt = opti.debug.value(X);
          U_opt = opti.debug.value(U);
          
          Uc(f,:) = U_opt;
          
          if strcmp(parms(i,f).type,'crossbridge')
              Xc(f,:) = X_opt(2,:);
          else
              Xc(f,:) = X_opt;
          end     
      end
      end

          
    % moved from 'find model inputs'
    function[U_smooth] = smoothen_solution(U,prob)

      U_smooth = nan(size(U));
      
      for f = 1:length(prob.freqs)
         
          % upsample
          tup = linspace(0,max(prob.t(f,:)),1000);
          Aup = interp1(prob.t(f,:), U(f,:), tup);
      
          fc = 30;
          fs = 1/mean(diff(tup));
          Wn = fc / (fs*.5);
          N = 2;
          [b,a] = butter(N, Wn,'low');
          
          Uup_smooth = filtfilt(b,a,Aup);
          
          % downsample
          U_smooth(f,:) = interp1(tup, Uup_smooth, prob.t(f,:));
      end
    end

    function[Xd] = dXfunc_shell(t, x, parms)
    
        r = parms.exp.A;
        Xd = cfxc.dXfunc(r, x, parms);
    end
    
    function[cost] = find_steadystate(X, parms)
%         parms.type = 'crossbridge';      
        
        [Xd] = cfxc.dXfunc(parms.exp.A, X, parms);
        cost = sum(Xd.^2);
    end
    
    function[cost] = find_CBrates(rates, parms, fv)
        
        parms.CB.f = rates(1);
        parms.CB.g(1) = rates(1);
        parms.CB.g(2:end) = rates(2:end);
        
        % CB force for Hill velocities
        FCB = cfxc.CB_force_velocity(fv.vHill, parms);
        
        % cost is sum-squared difference
        cost = sum((fv.FHill-FCB).^2);
    end
    
    
    function[n,p,q] = n_func(Q0, Q1, Q2, parms)

      if nargin < 4
          parms.xi = -3:.01:3;
      end
      
          eps = 1e-6;
      
          Q0c = max(Q0, eps); % correct because can't be zero

          p = Q1/Q0c; % Eq. 52
          q = sqrt(max(Q2/Q0c - (Q1/Q0c)^2, eps));  % Eq. 52
      
          % n
          n = Q0 ./ (sqrt(2*pi)*q) * exp(-((parms.xi-p).^2) / (2*q^2));  % Eq. 52, modified
    end

    function error = eulerIntegrator(x,xplus,uplus,dt)
      error = (xplus - x) - uplus*dt;
    end
    
    function[I] = TimTrapz(x, y)
      % I is the integral of y over x, following trapezoid rule
      dx = mean(diff(x));
      
      I = dx * (sum(y(1:end-1) + y(2:end))/2);
      
    end

    %%% all functions below were from main originally.
    %%% removed to allow testing/calling from other files.

    function[dX] = sim_muscle(t, x, parms)
        
        if parms.set.sim_mtc
            parms.exp.lmtc = x(end);
            xm = x(1:end-1); % muscle state
        else
            xm = x;
        end

        dX = xm(:);        
%       disp(t)
      Ca = x(1);

      if strcmp(parms.exp.stim_type,'u_func')
        u = parms.exp.u_func(t, parms);
      elseif strcmp(parms.exp.stim_type,'interp')
        u = interp1(parms.exp.t, parms.exp.U, t,'spline');
      elseif strcmp(parms.exp.stim_type,'max')
          u = 1;
      elseif strcmp(parms.exp.stim_type,'constant')
          u = parms.exp.A;
      end

      u(u<parms.ce.amin) = parms.ce.amin;

      % calcium dynamics
      dX(1,1) = parms.func.act(u, Ca, parms);
      %     Ca(Ca<parms.ce.amin) = parms.ce.amin; % bounds on the minimal activation

% if parms.isokinetic
%       parms.exp.lmtc = parms.Lfunc(t, parms);
%       parms.exp.vmtc = parms.Vfunc(t, parms);
% end

      if strcmp(parms.exp.stim_type,'max')
          dX(1,1) = 0;
      end
           
      if strcmp(parms.type, 'crossbridge')
        dX(2:end,1) = cfxc.crossbridge(Ca, xm(2:end), parms);
        
      elseif strcmp(parms.type, 'crossbridge_v2')
        dX(2:4,1) = cfxc.crossbridge_v2(Ca, xm(2:4), parms);
      
      elseif strcmp(parms.type, 'crossbridge_new')
        dX(2:end,1) = cfxc.crossbridge_new(Ca, xm(2:end), parms);
      
      elseif strcmp(parms.type, 'Huxley')
        dX(2:(parms.CB.nbins+2),1) = cfxc.crossbridge_original(Ca, xm(2:end), parms);
        
      elseif strcmp(parms.type, 'Hill-type')
        dX(2,1) = cfxc.Hill_type(Ca, xm(2), parms);

      elseif strcmp(parms.type, 'CaFaXC')
        Act = parms.func.sigmoid(Ca, parms);
%         Act = Ca;

        dX(2,1) = parms.func.fac(Act, xm(2), parms); % facilitation dynamics
        dX(3:end,1) = cfxc.crossbridge(xm(2), xm(3:end), parms); % cross-bridge dynamics
        
      elseif strcmp(parms.type, 'CaFaXC_v2')
        dX(2,1) = parms.func.fac(Ca, xm(2), parms); % facilitation dynamics
        dX(3:end,1) = cfxc.crossbridge_v2(xm(2), xm(3:end), parms); % cross-bridge dynamics        
      
      elseif strcmp(parms.type, 'CaFaC')
        dX(2,1) = parms.func.fac(Ca, xm(2), parms); % facilitation dynamics
        r = parms.func.sigmoid(Ca, parms);
        dX(3,1) = cfxc.Hill_type(r, xm(3), parms); % Hill-type dynamics      
      end
      
      % if we're simulating the lmtc
      if parms.set.sim_mtc
          dX(end+1,1) = parms.exp.vmtc;
      end
    end
    
    function[y,X] = get_model_output(t, x, parms)
        
         dX = nan(size(x))';
         y.D = nan(1,length(x));

         for i = 1:size(x,1)
            dX(:,i) = cfxc.sim_muscle(t(i), x(i,:), parms);
        end
        y.vce = dX(end,:)';

        y.lmtc = parms.exp.lmtc;
        y.vmtc = parms.exp.vmtc;
        
%         if parms.isokinetic == 1
%             y.lmtc = parms.Lfunc(t, parms);
%             y.mtmc = parms.Vfunc(t, parms);
%         end

        % states
        y.Ca = x(:,1);
        y.lce = x(:,end);
        y.Fpe = parms.func.fpe(y.lce, parms);
        y.lse = y.lmtc - y.lce;
        y.Fse = parms.func.fse((y.lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;
        y.Fce2 = y.Fse - y.Fpe;

        X = [y.Ca y.Fse y.lce];

        % contractile element force
        if strcmp(parms.type,'crossbridge') || strcmp(parms.type,'crossbridge_new')
            y.Fce = x(:,3) /(parms.CB.Xmax(2)/parms.ce.Fmax);
        elseif  strcmp(parms.type,'CaFaXC')
            y.Fce = x(:,4) /(parms.CB.Xmax(2)/parms.ce.Fmax);
            X(:,4) = x(:,2);
            y.Fa = x(:,2); 

         for i = 1:size(x,1)
            n = cfxc.n_func(x(i,3), x(i,4), x(i,5), parms.CB);
            n(parms.CB.xi > 1) = 0;
            y.D(i) = trapz(parms.CB.xi, n.*parms.CB.g_func(parms));
        end

        elseif  strcmp(parms.type,'Hill-type')
            y.Fce = y.Fse - y.Fpe;
        elseif  strcmp(parms.type,'CaFaC')
            y.Fce = y.Fse - y.Fpe;
            X(:,4) = x(:,2);
            y.Fa = x(:,2); 
        elseif strcmp(parms.type,'Huxley')
            gamma = parms.CB.h / (0.5 * parms.CB.s); % crossbridge to half-sarcomere
            alpha = 1 / (gamma * parms.ce.lceopt);
            dX = (y.lce - parms.exp.l0) * alpha;
            n = x(:,2:end-1);
            fmax_Huxley = parms.CB.f / (2*(parms.CB.f + parms.CB.g(1)));

            y.Fce = nan(size(Fse));
            for i = 1:size(x,1)
                xi = parms.CB.xi0(:) + dX(i);
                y.Fce(i) = trapz(xi(:), xi(:) .* n(i,:)') / fmax_Huxley * parms.ce.Fmax;
            end
        end
    end
    
function[Xd] = dXfunc(U, X, parms)
      U = U .* (U > 0);
%       X = X .* (X > 0);
      parms.set.fixed_velocity = 1;
          
      if strcmp(parms.type, 'crossbridge') || strcmp(parms.type, 'CaFaXC')
          % in this function, we treat velocity as input
          Xd = cfxc.crossbridge(U, X, parms);
      
      elseif strcmp(parms.type, 'crossbridge_new')
            Xd = cfxc.crossbridge_new(U, X, parms);
            
      else % simple first-order
           Xd = (U-X) / parms.ce.tau;
      end
end
    
    % moved from find_model_inputs
    function[beta,phi] = beta_phi_func(Q0,Q1,Q2,parms)

      if parms.CB.analytical
          
          if strcmp(parms.CB.ratefunc_type, 'Zahalak1981')
              beta = [parms.CB.f/2; parms.CB.f/3; parms.CB.f/4];
              phi = cfxc.analytical_solution_zahalak(Q0, Q1, Q2, parms);
          elseif strcmp(parms.CB.ratefunc_type, 'vanderZee2024')
              beta = 2./[1 2 3] .* parms.CB.f .* ((parms.CB.ps+parms.CB.w).^[1 2 3] - (parms.CB.ps-parms.CB.w).^[1 2 3]);
              phi = cfxc.analytical_solution_vanderzee(Q0, Q1, Q2, parms);
          end
      
      else
          % get distribution
          n = cfxc.n_func(Q0, Q1, Q2, parms.CB);
      
          beta0 = cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^0 .* parms.CB.f_func(parms.CB));   % Eq. 48
          beta1 = cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^1 .* parms.CB.f_func(parms.CB));   % Eq. 48
          beta2 = cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^2 .* parms.CB.f_func(parms.CB));   % Eq. 48
          

          if strcmp(parms.CB.ratefunc_type, 'Zahalak1981')
              phi0   = cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^0 .* parms.CB.f_func(parms.CB) .* n) + cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^0 .* parms.CB.g_func(parms.CB) .* n);   % Eq. 49+50
              phi1   = cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^1 .* parms.CB.f_func(parms.CB) .* n) + cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^1 .* parms.CB.g_func(parms.CB) .* n);   % Eq. 49+50
              phi2   = cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^2 .* parms.CB.f_func(parms.CB) .* n) + cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^2 .* parms.CB.g_func(parms.CB) .* n);   % Eq. 49+50
          else
              N = trapz(parms.CB.xi, n);
              phi0   = cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^0 .* parms.CB.f_func(parms.CB) .* N) + cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^0 .* parms.CB.g_func(parms.CB) .* n);   % Eq. 49+50
              phi1   = cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^1 .* parms.CB.f_func(parms.CB) .* N) + cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^1 .* parms.CB.g_func(parms.CB) .* n);   % Eq. 49+50
              phi2   = cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^2 .* parms.CB.f_func(parms.CB) .* N) + cfxc.TimTrapz(parms.CB.xi, parms.CB.xi.^2 .* parms.CB.g_func(parms.CB) .* n);   % Eq. 49+50
          end
          
          beta = [beta0 beta1 beta2];
          phi = [phi0 phi1 phi2];
      end
    end
    
    % moved from find_model_inputs
    function[beta, phi1, phi2, Rd] = beta_phi_func_v2(Q0,Q1,Q2,R,parms)
        % idea: could consider removing this and just adding additional
        % terms to the phi calculated by beta_phi_func (and adding Rd).
        % Would need to add g3 calculations, and set g3 = 0 if ripping.
        % Might be same with SRX state, need to discount crossbridges in
        % this state in the phi term.
                
        % for a model with a ripped state

        Q = [Q0 Q1 Q2];
        
        eps = 1e-6;

        p = Q(2)/max(Q(1), eps); % Eq. 52
        q = sqrt(max(Q(3)/max(Q(1), eps) - (Q(2)/max(Q(1), eps))^2, eps));  % Eq. 52

        c = [Q(1) ./ (sqrt(2*pi)*q) p 2*q^2];
       
        
        phi1 = nan(1,3);
        phi2 = nan(1,3);
        beta = nan(1,3);
    
      if parms.CB.analytical
          
            for i = 1:3
                % points were integrals is evaluated
                xi = [-parms.CB.K 0 parms.CB.dLcrit parms.CB.K];

                % most expensive line (50% CPT)
                IGs = parms.CB.gaussian.IG{i}(xi, c);

                if i == 1
                    Rd = (parms.CB.k * (IGs(4) - IGs(3)) -  (parms.CB.b * R) * ((parms.CB.ps+parms.CB.w) - (parms.CB.ps-parms.CB.w))); 
                end
                
                % beta: whatever is activation-dependent
                beta(i) = 1/i * (2 * parms.CB.f * 1) * ((parms.CB.ps+parms.CB.w)^i - (parms.CB.ps-parms.CB.w)^i);
                
                % phi1: whatever is activation-independent and SRX-independent          
                phi2(i) = parms.CB.g(1) * (IGs(3) - IGs(2)) ...
                        + parms.CB.g(2) * (IGs(2) - IGs(1)) ... % note: changed IGs(1) to IGs(2) in first half
                        + parms.CB.k(1) * (IGs(4) - IGs(3)) ...
                        - 1/i * (parms.CB.b * R) * ((parms.CB.ps+parms.CB.w)^i - (parms.CB.ps-parms.CB.w)^i);

                % phi1: whatever is activation-independent but SRX-dependent
                phi1(i) = 1/i * (2 * parms.CB.f * Q0) * ((parms.CB.ps+parms.CB.w)^i - (parms.CB.ps-parms.CB.w)^i);
            end
        
      end
    end

      function[beta, phi1, phi2, Rd] = beta_phi_func_v3(Q0,Q1,Q2,R,parms)

        % put bounds on the moments
        eps = [1e-6 1e-9];
        sigma_min = sqrt(eps(2)/eps(1)); % minimal s.d.

        Q0 = max(Q0, eps(1));
        Q1 = max(Q1, eps(1));
        Q2 = max(Q2, Q1^2/Q0+eps(2));

        % define some variables
        p = Q1/Q0; % mean of the distribution
        q = max(sqrt(max(Q2/Q0 - p^2,0)), sigma_min); % the standard deviation of the distribution
        
        % area, mean and scaled scared s.d.
        c1 = [Q0 p 2*q.^2];
        c2 = [1 parms.CB.ps 2*parms.CB.w.^2];
        Q = [Q0 Q1 Q2];

        % points were integrals is evaluated
        xi = [-parms.CB.K parms.CB.m+parms.CB.ps parms.CB.dLcrit(1) parms.CB.K];

        % analytical
        phi1 = nan(1,3);
        phi2s = nan(2,3);
        beta = nan(1,3);
        
%         parms.CB.K = 30;
%         x = linspace(-parms.CB.K, parms.CB.K, 1000);

        for i = 1:3

%             try
%             IGs_d = parms.CB.gaussian.IG{i}(xi, c1); % IGs for detachment 
            IGs_d = parms.CB.gaussian.IG2{i}(xi, c1, Q); % IGs for detachment 
%             catch
%             keyboard
%             end
%             IGs_a = parms.CB.gaussian.IG{i}(xi, c2); % IGs for attachment

            % breaking
            phi2s(1,i) = parms.CB.g1 * Q(i) ... % full integral
                        +parms.CB.g2 * (IGs_d(2) - IGs_d(1)) ...
                        +parms.CB.g3 * (IGs_d(4) - IGs_d(3));

            % ripping
%             phi2s(2,i) = +parms.CB.k  * (IGs_d(4) - IGs_d(3)) ...
%                          -parms.CB.b  * (IGs_a(4) - IGs_a(1)) * R(1);
            phi2s(2,i) = 0;

            % attaching
%             beta(1,i) =  parms.CB.f * (IGs_a(end) - IGs_a(1));
%             phi1(1,i) =  parms.CB.f * (IGs_a(end) - IGs_a(1)) * Q0;

            beta(1,i) =  parms.CB.f * parms.CB.gaussian.IGf{i}(c2);
            phi1(1,i) =  parms.CB.f * parms.CB.gaussian.IGf{i}(c2) * Q0;

%             betan(1,i) = trapz(x, x.^(i-1).*parms.CB.gaussian.G(x, c2) * parms.CB.f);
%             phi2n(1,i) = trapz(x, x.^(i-1).*parms.CB.gaussian.G(x, c1) * parms.CB.g1);
        end

        Rd = phi2s(2,1);
        phi2 = sum(phi2s);

      end


      function[beta, phi, Rd] = beta_phi_func_v4(Q0,Q1,Q2,R,SRX,parms)
        % idea: could consider removing this and just adding additional
        % terms to the phi calculated by beta_phi_func (and adding Rd).
        % Would need to add g3 calculations, and set g3 = 0 if ripping.
        % Might be same with SRX state, need to discount crossbridges in
        % this state in the phi term.
                
        % for a model with a ripped state

        Q = [Q0 Q1 Q2];
        
        eps = 1e-6;

        p = Q(2)/max(Q(1), eps); % Eq. 52
        q = sqrt(max(Q(3)/max(Q(1), eps) - (Q(2)/max(Q(1), eps))^2, eps));  % Eq. 52

        c = [Q(1) ./ (sqrt(2*pi)*q) p 2*q^2];
        
%         D = 1 - Q0 - R;
        
        phi = nan(1,3);
        beta = nan(1,3);
    
      if parms.CB.analytical
          
            for i = 1:3
                % points were integrals is evaluated
                xi = [-parms.CB.K 0 parms.CB.dLcrit parms.CB.K];

                if i == 1
                    % most expensive line (50% CPT)
                    IGs = parms.CB.gaussian.IG{1}(xi, c);
                    Rd = (parms.CB.k * (IGs(4) - IGs(3)) -  (parms.CB.b * R) * ((parms.CB.ps+parms.CB.w) - (parms.CB.ps-parms.CB.w))); 
                else
                    IGs = parms.CB.gaussian.IG{i}(xi, c, parms.CB.gaussian.G, parms.CB.gaussian.IG{1});
                end
                
                % beta: whatever is activation-dependent
                beta(i) = 1/i * (2 * parms.CB.f * 1) * ((parms.CB.ps+parms.CB.w)^i - (parms.CB.ps-parms.CB.w)^i);
                
                % phi: whatever is activation-independent          
                phi(i) = parms.CB.g(1) * (IGs(3) - IGs(2)) ...
                       + parms.CB.g(2) * (IGs(2) - IGs(1)) ... % note: changed IGs(1) to IGs(2) in first half
                       + parms.CB.k    * (IGs(4) - IGs(3)) ...
                       - 1/i * (parms.CB.b * R) * ((parms.CB.ps+parms.CB.w)^i - (parms.CB.ps-parms.CB.w)^i) ...
                       + 1/i * (2 * parms.CB.f * (R+Q0+SRX)) * ((parms.CB.ps+parms.CB.w)^i - (parms.CB.ps-parms.CB.w)^i);
            end
        
      end
    end
     
    function[beta, phi1, phi2, Rd] = beta_phi_func_v5(Q0,Q1,Q2,R,parms)

        % put bounds on the moments
        Q = [1e-6 1e-9];
        sigma_min = sqrt(eps(2)/eps(1)); % minimal s.d.

%         Q0 = max(Q0, eps(1));
%         Q1 = max(Q1, eps(2));
%         Q2 = max(Q2, Q1^2/Q0+eps(2));

        % define some variables
        p = Q1/Q0; % mean of the distribution
%         q = max(sqrt(max(Q2/Q0 - p^2,0)), sigma_min); % the standard deviation of the distribution
        q = sqrt(Q2/Q0 - p^2); % the standard deviation of the distribution

        % area, mean and scaled scared s.d.
        c1 = [Q0 p 2*q.^2];
        c2 = [1 0 2*parms.CB.w.^2];
        Q = [Q0 Q1 Q2];

        % points were integrals is evaluated
%         xi1 = [0 inf];
        k1 = [parms.CB.k11 parms.CB.k12];
        k2 = [parms.CB.k21 -parms.CB.k22];
%         k3 = [parms.k31 -parms.k32];

        % analytical
        phi1 = nan(1,3);
        phi2s = zeros(4,3);
        beta = nan(1,3);

        for i = 1:3

            % breaking
            phi2s(1,i) = parms.CB.gaussian.IGef{i}(c1,k1);
            phi2s(2,i) = parms.CB.g1 * Q(i); ... % full integral
            phi2s(3,i) = parms.CB.gaussian.IGef{i}(c1,k2);
%             phi2s(3,i) = parms.CB.g2 * parms.CB.gaussian.IGh{i+1}(c1); 

            % ripping
%             phir = parms.CB.gaussian.IGef{i}(c1,k2);
%             phir = 0;

%             if Q(1) <= 0
%                 phir = 0;
%             end

%             phi2s(4,i) = phir - parms.CB.b * parms.CB.gaussian.IGf{i}(c2) * R(1);
%             phi2s(4,i) = 0;
%             phi2s(2,i) = parms.CB.k * parms.CB.gaussian.IG{i}(inf,c1,k2) - parms.CB.b * parms.CB.gaussian.IGf{i}(c2) * R(1);

            % attaching
            beta(1,i) = parms.CB.f * parms.CB.gaussian.IGf{i}(c2);
            phi1(1,i) = parms.CB.f * parms.CB.gaussian.IGf{i}(c2) * Q(1);
        end

        Rd = phi2s(end,1);
        phi2 = sum(phi2s);

    end

    function[phi] = analytical_solution_zahalak(Q0, Q1, Q2,parms)          
        eps = 1e-6; 
    %     eps = 0;
        Q0c = max(Q0, eps);
        
        p = Q1/Q0c; % Eq. 52
        q = sqrt(max(Q2/Q0c - (Q1/Q0c)^2, eps));  % Eq. 52
        
        f1 = parms.CB.f(1);
        g1 = parms.CB.g(1);
        g2 = parms.CB.g(2);
        g3 = parms.CB.g(3);    
        
        % Jer
        jerf = @(tau) 1/2 * (1 + erf(tau/sqrt(2))); % translation of Zahalak's erf to the actual erf: see Sure hope this is right. % jer had *ci.sqrt(2) in the prefixing coefficient denominator.

        % following the J's from Zahalak 1981, A11.  
        J0 = @(tau,p,q) jerf(tau);
        J1 = @(tau,p,q) p^1*jerf(tau) -         q*exp(-tau^2/2) / sqrt(2*pi);
        J2 = @(tau,p,q) p^2*jerf(tau) - 2*p^1*  q*exp(-tau^2/2) / sqrt(2*pi) ...
            + q^2        * (jerf(tau) - tau*exp(-tau^2/2) / sqrt(2*pi));
        J3 = @(tau,p,q) p^3*jerf(tau) - 3*p^2*  q*exp(-tau^2/2) / sqrt(2*pi) ...
            + 3*p*q^2    * (jerf(tau) - tau*exp(-tau^2/2) / sqrt(2*pi)) ...
            - q^3*(2+tau^2) * exp(-tau^2/2) / sqrt(2*pi);
    
        % A-12. 
        phi0 = Q0*(g2*J0(-p/q,p,q) + (f1+g1)*( J1((1-p)/q,p,q) - J1(-p/q,p,q) )...
            + g1*(p             - J1((1-p)/q,p,q) ) + g3*(p                 - J1((1-p)/q,p,q) - 1 + J0((1-p)/q,p,q) ));
        phi1 = Q0*(g2*J1(-p/q,p,q) + (f1+g1)*( J2((1-p)/q,p,q) - J2(-p/q,p,q) )...
            + g1*(p^2+q^2       - J2((1-p)/q,p,q) ) + g3*(p^2+q^2           - J2((1-p)/q,p,q) - p + J1((1-p)/q,p,q)));
        phi2 = Q0*(g2*J2(-p/q,p,q) + (f1+g1)*( J3((1-p)/q,p,q) - J3(-p/q,p,q) )...
            + g1*(p^3+3*p*q^2   - J3((1-p)/q,p,q) ) + g3*(p^3 + 3*p*q^2     - J3((1-p)/q,p,q) - (p^2+q^2) + J2((1-p)/q,p,q)));
    
        phi = [phi0; phi1; phi2];

      end
      
      function[phi] = analytical_solution_vanderzee(Q0, Q1, Q2,parms) 
        Q = [Q0 Q1 Q2];

        eps = 1e-6;

        p = Q(2)/max(Q(1), eps); % Eq. 52
        q = sqrt(max(Q(3)/max(Q(1), eps) - (Q(2)/max(Q(1), eps))^2, eps));  % Eq. 52

        c = [Q(1) ./ (sqrt(2*pi)*q) p 2*q^2];


        phi = nan(1,3);
        IGs = nan(4,4);

        for i = 1:4
            % points were integrals is evaluated
            xi = [-parms.CB.K 0 parms.CB.dLcrit parms.CB.K];

            if i == 1
                % most expensive line (50% CPT)
                IGs(i,:) = parms.CB.gaussian.IG{1}(xi, c);
            else
                IGs(i,:) = parms.CB.gaussian.IG{i}(xi, c, parms.CB.gaussian.G, parms.CB.gaussian.IG{1});
            end
        end
        
        for i = 1:3
            % phi: whatever is activation-independent (note: this is a bit like Zahalak, 1986)                   
            phi(i)  = parms.CB.g(1)                     * (IGs(i,3) - IGs(i,2)) ... % g1 with x = 0-Lcrit
                    + (parms.CB.g(3)+parms.CB.g(1))     * (IGs(i+1,4) - IGs(i+1,3)) ... % g1+g3 for x > Lcrit
                    + parms.CB.g(2)                     * (IGs(i,2) - IGs(i,1)) ...; % g2 for x < 0
                    + 1/i * (2 * parms.CB.f * Q0) * ((parms.CB.ps+parms.CB.w)^i - (parms.CB.ps-parms.CB.w)^i); % f within binding region
        end
        
      end

    function[lce] = calc_length_from_force(F, parms)
        
        % method 1: ignore PE
        dlse = parms.func.lse(F/parms.ce.Fmax, parms) * parms.see.lse0;
        lcex = parms.exp.lmtc - dlse - parms.see.lse0; % lce along mtc
        
        % method 2: include PE, only works for certain functions
        
        % method 3: root-finding, is costly (not recommended)
        
        if parms.ce.pennation % assume constant thickness
            lce = sqrt(lcex.^2 + parms.ce.thickness^2);
        else
            lce = lcex; 
        end
    end
    
      function[Xd] = crossbridge_new(r, x, parms)

        %% retrieve states
        Q0 = x(1); 
        Q1 = x(2); 
        Q2 = x(3); 
        
        % we use Q1 here instead of Fcerel, otherwise difficult to
        % determine Xmax
%         DRX = r * (parms.CB.k1 + parms.CB.kF*Q1) / (parms.CB.k1 + parms.CB.kF*Q1 + parms.CB.k2);  
%         DRX = r;
        
        if length(x)<4
            R = 0;
        else
            R = x(4);
        end

        % rate constants
        parms.CB.f = parms.CB.scale_rates(r,parms.CB.f, parms.CB.mu);
        parms.CB.g = parms.CB.scale_rates(r,parms.CB.g, parms.CB.mu);

        % calculate beta and phi
%         [beta,phi1,phi2,Rd] = cfxc.beta_phi_func_v2(Q0,Q1,Q2,R,parms);
        [beta,phi1,phi2,Rd] = cfxc.beta_phi_func_v3(Q0,Q1,Q2,R,parms);
        %% determine length and filament overlap
        if ~parms.set.fixed_velocity
            % length scaling parameter
            gamma = parms.CB.h / (0.5 * parms.CB.s); % crossbridge to half-sarcomere
            alpha = 1 / (gamma * parms.ce.lceopt); % crossbridge to whole-muscle

            % either length is a quasi-state, or we need to calculate it
            if (length(x) > 5 && ~parms.thin_filament_coop) || (length(x) > 6)
                lce = x(end);
            else
                lce = cfxc.calc_length_from_force(Q1/(parms.CB.Xmax(2)/parms.ce.Fmax), parms);
            end

            if parms.set.optimum
                  a = parms.exp.a; % in case you assume optimum length
            else, a = parms.func.fce(lce, parms); % in case you use force-length relation
            end
       
       else
            a = parms.exp.a;
        end

        %% cooperative dynamics
        Ntot = a * r;

        if parms.thin_filament_coop
            Non = x(5);
        else
            Non = a*r;
        end

        % thin filament activation        
        if parms.thin_filament_coop
            
            % Campbell 2018
            Jon = 1/2 * parms.CB.kon  * (Ntot-Non)  * (1 + parms.CB.kcoop * Non/Ntot);
            Joff =      parms.CB.koff * (Non-Q0)    * (1 + parms.CB.kcoop * (Ntot-Non)/Ntot);

            Nond = Jon - Joff;
        
        else
            Nond = [];
        end

        if length(x)<6
            DRX = 1;
            SRXd = [];
        else
            SRX = x(6);
            DRX = 1 - SRX; % proportion of available crossbridges
            SRXd = parms.CB.k2 .* DRX - (parms.CB.k1 + parms.CB.kF*max(Q1,0)/r) .* SRX;  % note: force must be positive
        end

       %% length dynamics  
       if ~parms.set.fixed_velocity

            % calculate velocity
            Q1dot = DRX * (Non * beta(2) - phi1(2)) - phi2(2);
            [ks,kp] = cfxc.get_kse_kpe(lce, parms);
            [vce, u] = cfxc.enforce_forcerate_constraint(Q0,ks,kp,Q1dot, parms.exp.vmtc, parms.CB.Xmax(2), parms);

            if parms.set.no_tendon
                u = alpha * parms.exp.vmtc;
                vcerel = u * gamma;
                vce = vcerel * parms.ce.lceopt;
            end
        else
            u = parms.exp.u;
        end

        %% state derivates
        Qd0 = DRX * (Non * beta(1) - phi1(1)) - phi2(1);
        Qd1 = DRX * (Non * beta(2) - phi1(2)) - phi2(2) + 1 * u * Q0;
        Qd2 = DRX * (Non * beta(3) - phi1(3)) - phi2(3) + 2 * u * Q1;

        %% total state derivative vector
        if length(x) > 6
            Xd = [Qd0; Qd1; Qd2; Rd; Nond; SRXd; vce];
        else
            Xd = [Qd0; Qd1; Qd2; Rd; Nond; SRXd];
        end
        
        Xd = Xd(1:length(x));

      end
    
    function[Qd] = crossbridge(r, x, parms)
    % Zahalak derivative function with 3 states (the Qs)
    % based on Zahalak & Ma (1990), J. Biomech. Eng.
    
    %% retrieve states
    Q0 = x(1); 
    Q1 = x(2); 
    Q2 = x(3); 

    % rate constants
    parms.CB.f = parms.CB.scale_rates(r,parms.CB.f, parms.CB.mu);
    parms.CB.g = parms.CB.scale_rates(r,parms.CB.g, parms.CB.mu);

    % calculate beta and phi. beta includes all activation-dependent terms, phi includes the activation-independent terms
    [beta,phi] = cfxc.beta_phi_func(Q0,Q1,Q2,parms);

    %% length dynamics  
    if ~parms.set.fixed_velocity
        % length scaling parameter
        gamma = parms.CB.h / (0.5 * parms.CB.s); % crossbridge to half-sarcomere
        alpha = 1 / (gamma * parms.ce.lceopt); % crossbridge to whole-muscle

    % either length is a quasi-state, or we need to calculate it
    if (length(x) > 3 && ~parms.thin_filament_coop) || length(x) > 4
        lce = x(end);
    else
        lce = cfxc.calc_length_from_force(Q1/(parms.CB.Xmax(2)/parms.ce.Fmax), parms);
    end

    if parms.set.optimum
        a = parms.exp.a; % in case you assume optimum length
    else, a = parms.func.fce(lce, parms); % in case you use force-length relation
    end

    else
     u = parms.exp.u;
     a = parms.exp.a;
    end

    if parms.thin_filament_coop
        Non = x(4);
    else
        Non = a*r;
    end

   if ~parms.set.fixed_velocity
        % calculate velocity
        Q1dot = Non * beta(2) - phi(2); % velocity-independent Q1dot
        [ks,kp] = cfxc.get_kse_kpe(lce, parms);
        [vce, u] = cfxc.enforce_forcerate_constraint(Q0,ks,kp,Q1dot, parms.exp.vmtc, parms.CB.Xmax(2), parms);

        if parms.set.no_tendon
            u = alpha * parms.exp.vmtc;
            vcerel = u * gamma;
            vce = vcerel * parms.ce.lceopt;
        end
    else
        vce = [];
    end

    if parms.thin_filament_coop
        % Campbell 2018 adjusted
        Ntot = a * r;
        Jon = 1/2 * parms.CB.kon  * (Ntot-Non)       * (1 + parms.CB.kcoop * Non);
        Joff =      parms.CB.koff * (Non-Q0)    * (1 + parms.CB.kcoop * (Ntot-Non));

        Nond = Jon - Joff;
    else
        Nond = [];
    end

    %% state derivates
    Qd0 = Non * beta(1) - phi(1);
    Qd1 = Non * beta(2) - phi(2) + 1 * u * Q0;
    Qd2 = Non * beta(3) - phi(3) + 2 * u * Q1;
    
    % total state derivative vector
    Qd = [Qd0; Qd1; Qd2; Nond; vce];
 
    end
    
    function[Qd] = crossbridge_v2(r, x, parms)
    % Zahalak derivative function with 3 states (the Qs)
    % based on Zahalak & Ma (1990), J. Biomech. Eng.
    % here we have length as a state

    %% retrieve states
    lce = x(1); 
    Q0 = x(2); 
    Q2 = x(3); 
    
    % determine force from length
    lse = parms.exp.lmtc - lce;
    dlse_rel = (lse - parms.see.lse0) / parms.see.lse0;
    Ft = parms.ce.Fmax * parms.func.fse(dlse_rel, parms);
    Fp = parms.func.fpe(lce, parms);
    Fm = Ft - Fp;
    Q1 = Fm * parms.CB.Xmax(2)/parms.ce.Fmax; % note: delta = fmax/Fmax (omgekeerd)
    
    % rate constants
    parms.CB.f = parms.CB.scale_rates(r,parms.CB.f, parms.CB.mu);
    parms.CB.g = parms.CB.scale_rates(r,parms.CB.g, parms.CB.mu);

    % calculate beta and phi
    [beta,phi] = cfxc.beta_phi_func(Q0,Q1,Q2,parms);
      
    % length scaling parameter
    gamma = parms.CB.h / (0.5 * parms.CB.s); % crossbridge to half-sarcomere

    %% length dynamics  
    if ~parms.set.fixed_velocity
    if parms.set.optimum
        a = parms.exp.a; % in case you assume optimum length
    else, a = parms.func.fce(lce, parms); % in case you use force-length relation
    end

    % calculate velocity
    Q1dot = a * r * beta(2) - phi(2); % velocity-independent Q1dot
   [ks,kp] = cfxc.get_kse_kpe(lce, parms);
    fmax = parms.CB.Xmax(2);
    [vce, u] = cfxc.enforce_forcerate_constraint(Q0,ks,kp,Q1dot, parms.exp.vmtc, fmax, parms);

    if parms.set.no_tendon
        alpha = 1 / (gamma * parms.ce.lceopt); % crossbridge to whole-muscle
        u = alpha * parms.exp.vmtc;
        vcerel = u * gamma;
        vce = vcerel * parms.ce.lceopt;        
    end

    else

     u = parms.exp.u;
     a = parms.exp.a;
    end

    %% state derivates
    Qd0 = a * r * beta(1) - phi(1);
    Qd2 = a * r * beta(3) - phi(3) + 2 * u * Q1;
    
    % total state derivative vector
    Qd = [vce; Qd0; Qd2];

    end
    
    function[Qd] = crossbridge_original(r, x, parms)

    % original Huxley (1957), similar to implemented by Lemaire et al. (2016)

    %% retrieve states
    n = x(1:end-1);
    lce = x(end);

    n = n(:);

    %% determine force-length
    if parms.set.optimum
    a = parms.exp.a; % in case you assume optimum length
    else, a = parms.func.fce(lce, parms); % in case you use force-length relation
    end

    %% recalc some parms
    fmax = parms.CB.f / (2*(parms.CB.f + parms.CB.g(1)));
    gamma = parms.CB.h / (0.5 * parms.CB.s); % crossbridge to half-sarcomere
    alpha = 1 / (gamma * parms.ce.lceopt);

    %% determine elastic stiffnesses
%     kp = parms.func.kpe(lce, parms); % [N/m]
%     dlse_rel = ((parms.exp.lmtc - lce)-parms.see.lse0) / parms.see.lse0;
%     ks = parms.func.kse(dlse_rel,parms) * (parms.ce.Fmax/parms.see.lse0); % [N/m]

    % rate constants
    parms.CB.f = parms.CB.scale_rates(r, parms.CB.f, parms.CB.mu); 
    parms.CB.g = parms.CB.scale_rates(r, parms.CB.g, parms.CB.mu); 

    % displacement from start
    dX = (lce - parms.exp.l0) * alpha;
    xi = parms.CB.xi0 + dX;
    iRel = ((xi(:) < 2) & (xi(:) > -1)) | (abs(n(:)) > 1e-16);

    % only select relevant portion
    parms.CB.xi = xi(iRel);
    n = n(iRel);

    beta = parms.CB.f_func(parms);
    phi = (parms.CB.f_func(parms) + parms.CB.g_func(parms)) .* n';

    ndot0 = a * r * beta - phi;

    Q0 = trapz(parms.CB.xi, n);
    Q1dot = trapz(parms.CB.xi, ndot0 .* parms.CB.xi);

   [ks,kp] = cfxc.get_kse_kpe(lce, parms);
    vce = cfxc.enforce_forcerate_constraint(Q0,ks,kp,Q1dot, parms.exp.vmtc, fmax, parms);
   

      if parms.set.no_tendon
          u = alpha * parms.exp.vmtc;
        vcerel = u * gamma;
        vce = vcerel * parms.ce.lceopt;          
      end

      % estimate spatial derivative
%       dndx = cfxc.grad5(n(:), mean(diff(parms.CB.xi)));
      
      ndot = zeros(size(parms.CB.xi0));

      % state derivates
      ndot(iRel) = a * r * beta(:) - phi(:);
      
      % dimensionless
%       vce = u / alpha;

      Qd = [ndot(:); vce(:)];

    end
    
    function[vce, u] = enforce_forcerate_constraint(Q0,ks,kp,Q1dot, vmtc, fmax, parms)
        % ks and kp should be relative to lceopt

        gamma = parms.CB.h / (0.5 * parms.CB.s); % crossbridge to half-sarcomere
       
        % calculate velocity that assures that forces are compatible at next time step
        vcerel = (vmtc/parms.ce.lceopt * ks - Q1dot/fmax) ...
        / (ks + kp + Q0/(fmax*gamma));

        % back to crossbridge
        u = vcerel / gamma;
        vce = vcerel * parms.ce.lceopt;
    end

    function[ks,kp] = get_kse_kpe(lce, parms)

        % determine elastic stiffnesses
        dlse_rel = ((parms.exp.lmtc - lce) - parms.see.lse0) / parms.see.lse0;
        kp = parms.func.kpe(lce, parms) * (parms.ce.lceopt / parms.ce.Fmax);
        ks = parms.func.kse(dlse_rel,parms) * (parms.ce.lceopt/parms.see.lse0); % expressed relative to lceopt

    end


    function[Vce] = Hill_type(a, lce, parms)

      Lse = parms.exp.lmtc - lce;
      Fse = parms.func.fse((Lse-parms.see.lse0)/parms.see.lse0,parms) * parms.ce.Fmax;
      
      a(a<parms.ce.amin) = parms.ce.amin;

      if Lse <= parms.see.lse0
        Fse = 0;
      end

      Fpe = parms.func.fpe(lce, parms);

      Fce = Fse - Fpe;
      Fce(Fce < 0) = 0;

      if parms.set.optimum 
          Fisom = parms.exp.a; % in case you assume optimum length
      else, Fisom = parms.func.fce(lce, parms); % in case you use force-length relation
      end

      % force-velocity
      Vce = parms.func.fv(a, Fce/parms.ce.Fmax, Fisom, parms) * parms.ce.lceopt;

      if isnan(Vce), keyboard
      end
    end

    function[parms] = calc_l0(parms)
    % lmtc
    parms.exp.lmtc = parms.func.lmtc(parms.exp.phi, parms);

    % CE length (assuming slack elastic tissues)
    lce0 = parms.exp.lmtc - parms.see.lse0;

    % assume minimal excitation
    parms.exp.A = parms.ce.amin;

    % if PE is exerting force, find length at which elastic forces are equal
    opt = optimset('TolFun',1e-20);
    if lce0 > (parms.pe.lpe0*parms.ce.lceopt)
    parms.exp.l0 = fminsearch(@(L0) cfxc.find_l0(L0,parms), lce0, opt);
    else % else, SE is also slack
    parms.exp.l0 = lce0;
    end

    % now we also know the myofilament overlap force-length parameter
    parms.exp.a = parms.func.fce(parms.exp.l0, parms);
  end
  
    function[parms] = calc_x0(parms)

   % calc l0
    parms = cfxc.calc_l0(parms);

    % now we can compute the crossbridge distribution for a certain (minimum) activation
    parms.exp.u = 0;

    % non-linear inequality constraint
    nlincon_ineq = @(X) (X(2)/X(1))^2 - X(3)/X(1); % < 0
    nlincon_func = @(X) deal(nlincon_ineq(X),[]); 
    LB = -inf * ones(size(parms.CB.Xmax));
    LB(1) = 1e-6; % Q0 can't be 0
    LB(4) = 1e-6; % R can't be <0
    LB(5) = 1e-6; % DRX can't be <0
    UB = ones(1,6);

    % isometric
    opt = optimset('Display','off');
    [X0, fval] = fmincon(@(X) cfxc.find_steadystate(X, parms), parms.ce.amin * parms.CB.Xmax, [],[],[],[],LB(1:length(parms.CB.Xmax)),UB(1:length(parms.CB.Xmax)), nlincon_func, opt);
    
    % initial conditions
    parms.exp.x0 = [parms.ce.amin X0 parms.exp.l0];
    
    if fval > 1e-3
        parms.type = 'crossbridge';
        parms.set.fixed_velocity = 1;
        parms.exp.stim_type = 'constant';
        parms.exp.A = parms.ce.amin;
        [t,x] = ode113(@cfxc.sim_muscle, [0 5], [parms.ce.amin parms.ce.amin * parms.CB.Xmax], parms.set.odeopt, parms);
        X0 = x(end,:);
        parms.exp.x0 = [X0 parms.exp.l0];
    end
    

    % calc Xmax
    opt = optimset('Display','on');
    parms.exp.A = 1;
    parms.exp.a = 1;
    [parms.CB.Xmax, fval] = fmincon(@(X) cfxc.find_steadystate(X, parms), parms.CB.Xmax, [],[],[],[],LB(1:length(parms.CB.Xmax)),[],nlincon_func, opt);
      
    % simulate if result is not satisfying
    if fval > 1e-3
        [t,x] = ode113(@cfxc.sim_muscle, [0 5], parms.exp.x0(1:end-1), parms.set.odeopt, parms);
        parms.CB.Xmax = x(end,2:end);
    end
    
    end

    function[cost] = find_l0(L0, parms)

      Fpe = parms.func.fpe(L0,parms);

      Fce =  parms.exp.A * parms.ce.Fmax * parms.func.fce(L0, parms);

      lse = parms.exp.lmtc - L0;
      Fse = parms.func.fse((lse-parms.see.lse0)/parms.see.lse0, parms) * parms.ce.Fmax;

      % Fse = interp1(parms.Lses, parms.Fses, parms.exp.lmtc-L0);

      cost = 100*((Fpe + Fce) - Fse).^2;
    end

    function yp=grad5(y,dx)
      % function yp=grad5(y,dx)
      % 031091 KvS
      % purpose: calculation of 5 point derivative
      % inputs : y : vector containing signal as a function of x
      %          dx: change in x between steps
      % output : yp: derivative of y with respect to x
       
      if nargin~=2
        disp('OOPS: grad5 expects 2 input arguments')
        disp('      hit any key to continue');pause
        return
      end
      
      %y=y(:);
      [nrow,ncol]=size(y);
      
      
      yp(1,1:ncol)=zeros(1,ncol);
      yp(2,1:ncol )=zeros(1,ncol);
      yp(nrow-1,1:ncol)=zeros(1,ncol);
      yp(nrow,1:ncol)=zeros(1,ncol);
      
      yp(1,1:ncol)=(y(2,:)-y(1,:))/dx;
      yp(2,1:ncol)=(y(3,:)-y(1,:))/(2*dx);
      yp(nrow-1,1:ncol)=(y(nrow,:)-y(nrow-2,:))/(2*dx);
      yp(nrow,1:ncol)=(y(nrow,:)-y(nrow-1,:))/dx;
      
      coef=[1 -8 0 8 -1];
      for i=3:nrow-2;
        yp(i,:)=(coef*y(i-2:i+2,:))/(12*dx);
      end;
      
      %if flip; y=y';yp=yp';end;
    end
    
    function[F,n,F2] = CB_force_velocity(v, parms)

        gamma = parms.CB.h / (0.5 * parms.CB.s); % crossbridge to half-sarcomere
        
        Fiso = 1;
        Fcon = (v < 0) .* ((1 - exp((parms.CB.f+parms.CB.g(1))*gamma./(2*v)))   .* (2*v ./ ((parms.CB.f+parms.CB.g(1))*gamma) - 2.*v.^2 / (gamma^2*parms.CB.g(2)^2)));
        Fecc = -(v > 0) .* ((1 - exp(-(parms.CB.f+parms.CB.g(1))*gamma./(2*v))) .* (2*v ./ ((parms.CB.f+parms.CB.g(1))*gamma) + 2*v ./ ((parms.CB.g(1)+parms.CB.g(3))*gamma) ...
                .* (-parms.CB.g(3)*sqrt(pi)*sqrt(gamma)./(sqrt((parms.CB.g(1)+parms.CB.g(3)).*2*v)) .*  exp((parms.CB.g(1)^2*gamma)./((parms.CB.g(1)+parms.CB.g(3)).*2.*v)) ...
                .* (1 - erf(parms.CB.g(1).*sqrt(gamma) ./ (sqrt((parms.CB.g(1)+parms.CB.g(3)).*2*abs(v))))) - 1)));
        
       Fcon(isnan(Fcon)) = 0;
       Fecc(isnan(Fecc)) = 0;
            
        F = Fiso + Fcon + Fecc;
            
        if nargout > 1
            % calculate n(x) (Eq. A9)
            x = parms.CB.xi * parms.CB.h;
            h = parms.CB.h;
            f1 = parms.CB.f(1)/parms.CB.h;
            g1 = parms.CB.g(1)/parms.CB.h;
            g3 = parms.CB.g(3)/parms.CB.h;
            g0 = parms.CB.g(2);

            v = (v * parms.CB.s) / 2;

            n = zeros(length(v), length(parms.CB.xi));
            F2 = nan(size(v));

            Fiso = f1 / (2*(f1+g1));

            for i = 1:length(v)
                if v(i) == 0
                    n(i,x>=0 & x < h) = f1 / (f1 + g1);

                elseif v(i) < 0

                    n(i,x>=0 & x < h) = f1 / (f1 + g1) * (1 - exp(-(g1+f1)*(x(x>=0 & x < h).^2-h.^2) / (2*v(i))));
                    n(i,x<0) = (f1 / (f1 + g1)) * (1 - exp((g1+f1)*h^2 / (2*v(i)))) * exp((-g0.*x(x<0)) / v(i));

                elseif v(i) > 0
                    n(i,x>=0 & x < h) = f1 / (f1 + g1) * (1 - exp(-(g1+f1)*(x(x>=0& x < h).^2) / (2*v(i))));


                    C = (exp(-(g1+f1)*h.^2 / (2*v(i)))-1) * -exp((g1-g3)*h^2 / (2*v(i)));
                    n(i,x>=h) = f1 / (f1 + g1) * C * exp((-((g1+g3)*x(x>=h).^2 - 2*g3*x(x>=h)*h)) / (2*v(i)));
                end

                   F2(i) = trapz(parms.CB.xi, parms.CB.xi .* n(i,:)) / Fiso;
            end
        end
    
 
    end

    % this function assumes CaFaXC 
    function dstate = sim_segment(t,state,parms)
       
      
      phi_deg       = state(end-1);
      omega_deg     = state(end); 
     
      phi_rad   = phi_deg * pi/180;
      omega_rad = omega_deg * pi/180;

      % larger angle -> muscle shortening
      dmtcdt    = -parms.mtc.r * omega_rad; %
      
      % before calling sim_muscle, we need to set two paramters:
      %parms.exp.vmtc
      %parms.exp.lmtc
      lmtc_cur    = parms.func.lmtc(-phi_deg,parms);
      parms.exp.vmtc  = dmtcdt;
      parms.exp.lmtc  = lmtc_cur; 

      % muscle derivative
      dmusdt    = cfxc.sim_muscle(t, state(1:end-2),parms);
      
      % state to output
      [y,X] = cfxc.get_model_output(t, state(1:(end-2))', parms);

      % pendulum dynamics
      r_gyr     = 0.2;
      M         = 5;
      I         = M*r_gyr^2; %ml^2
 
      g         = 9.81;
      r_g       = r_gyr*cos(phi_rad);
      tau_g     = -r_g * M*g;
      tau_mus   = y.Fse*parms.mtc.r;
      
      % passive joint torque
%       tau_p = -(100 * joint_rad.^2) .* (joint_rad < 0);
      
      alpha_rad = (tau_mus + tau_g)/I;
      alpha_deg = alpha_rad * 180/pi;

      if strcmp(parms.mode,'pre-movement')
          alpha_deg = -parms.A * (2*pi*parms.f)^2*cos(2*pi*t*parms.f);
          omega_deg = -parms.A * 2*pi*parms.f*sin(2*pi*t*parms.f);
      end
          
          
%       ddjointdt2= -ddworlddt2;
      dstate    = [dmusdt;omega_deg;alpha_deg];
    end

    function clcs = demo_sim_leg(parms)
      
      % knee extension movement, begin knee bent straight down.
      phi_rad_0 = 0;

      if isempty(parms)
      load('quad_parms.mat','parms');
      end
      
      parms.exp.phi = phi_rad_0 * 180/pi;
      
      parms = cfxc.calc_x0(parms); 

      parms.exp.stim_type = 'constant';
      parms.exp.A = .02;

      % parms has too many parameters; here we specify explicity for some code control.
      parms.type = 'CaFaXC';
      parms.set.optimum = 0; % a flag throughout the code to assume lce=lceopt.
      parms.set.sim_mtc = 0; % 0->compute the velocity, rather than commanding it.
        
      % initial state
      X0 = [parms.exp.x0(1), parms.exp.x0(1:end),phi_rad_0,0];
      
      [t_s,state_s] = ode113(@cfxc.sim_segment, [0 10], X0, parms.set.odeopt, parms);
      clcs = cfxc.calcs_facilitation_sim(t_s,state_s,parms);

    end

    function out = calcs_facilitation_sim(t_s,state_s,parms)
      out.ca            = state_s(:,1);
      out.act           = state_s(:,2);
      out.Q0            = state_s(:,3);
      out.Q1            = state_s(:,4);
      out.Q2            = state_s(:,5);
      out.lce           = state_s(:,6);
      out.phi_joint_rad = state_s(:,7);
      
      out.lmtc        = parms.func.lmtc(out.phi_joint_rad*180/pi,parms);
      out.lse         = out.lmtc-out.lce;
      out.lse_frac_ex = (out.lse-parms.see.lse0)/parms.see.lse0;
      out.F_se        = parms.func.fse(out.lse_frac_ex,parms)*parms.ce.Fmax;
      out.F_pe        = parms.func.fpe(out.lce,parms);
      out.F_mus       = out.Q1/(parms.CB.Xmax(2)/parms.ce.Fmax);
      out.F_diff      = out.F_se-out.F_mus-out.F_pe;

      figure();
      subplot(221)
      plot(t_s, out.ca); hold on
      plot(t_s, out.act); hold on
      xlabel('time (s)');
      ylabel('activation')
      
      subplot(222)
      plot(t_s, 180-out.phi_joint_rad*180/pi); hold on
      xlabel('time (s)');
      ylabel('joint angle (deg)')
      ylim([20 180])
      
      subplot(223);
      plot(t_s,out.F_mus);hold on;
      plot(t_s,out.F_se);
      plot(t_s,out.F_pe);
      plot(t_s,out.F_pe+out.F_mus,'--');
      ylabel('force (N)')
      xlabel('time (s)')
      legend('CE','SE','PE','CE+PE','location','best');
      legend boxoff
            
      subplot(224);
      plot(t_s,out.lce);hold on;
      plot(t_s,out.lse);
      plot(t_s,out.lmtc);
      plot(t_s,repmat(parms.see.lse0,length(t_s),1),'k--');
      plot(t_s,repmat(parms.pe.lpe0*parms.ce.lceopt,length(t_s),1),'k:');
      ylabel('lengths (m)')
      xlabel('time (s)');
      
      legend('CE','SE','MTC','SE slack','PE slack','location','best');
      legend boxoff
      
      for i = 1:4
          subplot(2,2,i)
%           axis tight
          box off
          xlim([0 max(t_s)])
      end

    end
    
    
        function[fl] = evaluate_force_length(phis, parms,show)
    Lces = linspace(0, 2*parms.ce.lceopt, 1000);
    Fces = linspace(0,.99*parms.ce.Fasymp*parms.ce.Fmax, 1000);

    fl.Fses = linspace(0, parms.ce.Fmax,1000);
    fl.Lses = parms.see.lse0 * parms.func.lse(fl.Fses/parms.ce.Fmax, parms) + parms.see.lse0;

    fl.Lpes = linspace(parms.pe.lpe0*parms.ce.lceopt, parms.pe.lpemax*parms.ce.lceopt, 1000);
    fl.Fpes = parms.func.fpe(fl.Lpes, parms);
        
    fl.phis = phis;
    fl.Fse = nan(length(fl.phis),2);
    fl.Fpe = nan(length(fl.phis),2);
    fl.Fces = nan(length(fl.phis),2);
    fl.Lces = nan(length(fl.phis),2);
    fl.lSE = nan(length(fl.phis),2);

    opt = optimset('display','iter', 'TolX',1e-12,'TolFun',1e-12);


    for p = 1:length(fl.phis)
        parms.exp.lmtc = parms.func.lmtc(fl.phis(p), parms);

        for i = 1:2

            if i == 1 % passive
%                 parms.exp.A = parms.ce.amin;
                parms.exp.A = 0;
                if (parms.exp.lmtc-parms.see.lse0) > (parms.pe.lpe0*parms.ce.lceopt)
                    fl.Lces(p,i) = fmincon(@(L0) cfxc.find_l0(L0,parms), parms.exp.lmtc - parms.see.lse0, 1, (parms.exp.lmtc-parms.see.lse0),...
                        [],[],[],[],[], opt);
                else
                    fl.Lces(p,i) = parms.exp.lmtc - parms.see.lse0;
                end

            elseif i == 2 % active
                parms.exp.A = 1;
                fl.Lces(p,i) = fmincon(@(L0) cfxc.find_l0(L0,parms), fl.Lces(p,1)-.01, 1, fl.Lces(p,1),...
                        [],[],[],[],[], opt);
            end

            fl.lSE(p,i) = parms.exp.lmtc - fl.Lces(p,i);    
            fl.Fse(p,i) = parms.ce.Fmax * parms.func.fse((fl.lSE(p,i)-parms.see.lse0)/parms.see.lse0, parms);
            fl.Fpe(p,i) = parms.func.fpe(fl.Lces(p,i), parms);
            fl.Fces(p,i) = parms.exp.A * parms.func.fce(fl.Lces(p,i), parms) * parms.ce.Fmax;
        end


    end

if show
    %% Figure 2: Muscle level
    if ishandle(2), close(2); end
    figure(2)
    subplot(221);
    plot(Lces, parms.func.fce(Lces, parms)*parms.ce.Fmax); 
    xline(parms.ce.lceopt,'k--')

    subplot(222);
    plot(parms.func.fv(1, Fces/parms.ce.Fmax, 1, parms), Fces);
    xlim(parms.ce.vmaxrel * [-1 1])
    yline(parms.ce.Fmax,'k--')
    xline(0,'k--')

    subplot(223);
    plot(fl.Lses, fl.Fses); hold on; box off

    yline(parms.see.sefm*parms.ce.Fmax,'k--')
    xline(parms.see.sexm*parms.see.lse0 + parms.see.lse0,'k--')

    subplot(224);
    plot(fl.Lpes, fl.Fpes);

    titles = {'CE force-length', 'CE force-velocity', 'SE force-length', 'PE force-length'};
    xlabels = {'l_M (m)', 'v_{CE} (m/s)', 'l_{SE} (m)', 'l_M (m)'};

     for i = 1:4
         subplot(2,2,i)
         ylabel('Force (N)')
         title(titles{i})
         xlabel(xlabels{i})
         box off
         hold on
     end
 
     %% Figure 2: Joint level
    if ishandle(3), close(3); end
    figure(3)

    subplot(221)
    plot(fl.phis,fl.Fse(:,1)/parms.ce.Fmax); hold on; box off
    plot(fl.phis,fl.Fces(:,1)/parms.ce.Fmax,'--'); hold on
    plot(fl.phis,fl.Fpe(:,1)/parms.ce.Fmax,'--'); hold on
    xlabel('Angle (deg)'); ylabel('Force (0-1)');
    legend('Fse','Fce','Fpe','location','best')
    legend boxoff
    title('Passive force-angle')
    axis([min(fl.phis) max(fl.phis) 0 1])

    subplot(222)
    plot(fl.phis,fl.Fse(:,2)/parms.ce.Fmax); hold on; box off
    plot(fl.phis,fl.Fces(:,2)/parms.ce.Fmax,'--'); hold on
    plot(fl.phis,fl.Fpe(:,2)/parms.ce.Fmax,'--'); hold on
    xlabel('Angle (deg)'); ylabel('Force (0-1)');
    legend('Fse','Fce','Fpe','location','best')
    legend boxoff
    title('Active force-angle')
    axis([min(fl.phis) max(fl.phis) 0 1])

    subplot(223)
    plot(fl.phis, fl.Lces); box off
    xlabel('Angle (deg)'); ylabel('Length (m)');
    legend('Passive', 'Active','location','best')
    legend box off
    title('Contractile element (CE) length')
    xlim([min(fl.phis) max(fl.phis)])

    subplot(224)
    plot(fl.phis, -diff(fl.Lces,[],2)); box off
    xlabel('Angle (deg)'); ylabel('\Delta Length (m)');
    title('CE length change')
    xlim([min(fl.phis) max(fl.phis)])
    
    colors = [linspace(0,1,length(fl.phis))', zeros(length(fl.phis),2)];
    symb = '.x';

    %% Figure 1
    figure(2)
    for p = 1:length(fl.phis)
        parms.exp.lmtc = parms.func.lmtc(fl.phis(p), parms);

        for i = 1:2

        subplot(223); hold on; box off
        plot(fl.lSE(p,i), fl.Fse(p,i),symb(i),'color',colors(p,:))

        subplot(224); hold on; box off
        plot(fl.Lces(p,i), fl.Fpe(p,i),symb(i),'color',colors(p,:))

        subplot(221); hold on; box off
        plot(fl.Lces(p,i), fl.Fces(p,i),symb(i),'color',colors(p,:))
        subplot(221);
        plot(fl.Lces(p,:), fl.Fces(p,:),'-','color',colors(p,:))

        end
    end

end
    %
    end
  
    function[parms, fv] = fit_CB_on_Hill(parms, fv, type)
       
        % settings and initial guess
        rates0 = [150 1000 100];
        
        if strcmp(type, 'fit_CB')
            MaxIter = 1e3;
        else
            MaxIter = 1e1;
        end
        
        fopt = optimset('display','iter', 'TolFun', 1e-1, 'TolX', 1e-1, 'MaxIter', MaxIter);
    
        if strcmp(type, 'fit_CB') || strcmp(type, 'fit_DM')
            % find the rates for original model
            CBrates = fminsearch(@(rates) cfxc.find_CBrates(rates, parms, fv), rates0, fopt);

            % rates
            parms.CB.f = CBrates(1);
            parms.CB.g = CBrates;
        
        end

        % Huxley-type
        fv.FCB = zeros(length(fv.vHill),2);
        [fv.FCB(:,1),fv.n,fv.FCB(:,2)] = cfxc.CB_force_velocity(fv.vHill, parms);
        
        if strcmp(type, 'fit_DM')
            % adjust rates to DM approximation (optional)
            DMrates = fminsearch(@(CBrates) cfxc.find_DMrates(CBrates, parms, fv), CBrates, fopt);

            % rates
            parms.CB.f = DMrates(1);
            parms.CB.g = DMrates;
        end
    
        % evaluate
        [fv.FCB(:,3), parms] = cfxc.evaluate_DM(fv.vHill, parms);
    end
    
    function[cost] = find_DMrates(rates, parms, fv)
    
        parms.CB.f = rates(1);
        parms.CB.g = rates;
              
        % only consider a portion
        ID = 1:2:length(fv.vHill);
        vHill = fv.vHill(ID);
        
        % get DM force
        FCB = cfxc.evaluate_DM(vHill, parms);
        
        % cost is sum-squared difference
        cost = sum((fv.FHill(ID)-FCB').^2);
        
    end
    
   
function[FCB, parms, kCB] = evaluate_DM(v, parms)

    % non-linear inequality constraint
    nlincon_ineq = @(X) (X(2)/X(1))^2 - X(3)/X(1); % < 0
    nlincon_func = @(X) deal(nlincon_ineq(X),[]); 
    LB = [1e-6 -inf -inf 1e-6 1e-6 1e-6]; % Q0 can't be 0

    % isometric conditions
    parms.exp.u = 0;
    parms.exp.a = 1; % at optimum length
    parms.exp.A = 1; % maximal excitation
    parms.exp.stim_type = 'max';
    parms.set.fixed_velocity = 1;
    
    % isometric
    opt = optimset('Display','off');
    X0 = parms.CB.Xmax;
    
    [parms.CB.Xmax, fval] = fmincon(@(X) cfxc.find_steadystate(X, parms), X0, [],[],[],[],LB(1:length(X0)),[],nlincon_func, opt);
    
    if fval > 1e-3
        disp('Could not solve isometric steady-state, simulating ...')
        [t,x] = ode113(@cfxc.sim_muscle, [0 5], [1 X0], parms.set.odeopt, parms);
        parms.CB.Xmax = x(end,2:end);
    end
    
    % concentric
    id = v<=0;
    vs = flip(v(id));

    Xmax = nan(length(vs), length(X0));
    cost = nan(length(vs), 1);

    for i = 1:length(vs)
        parms.exp.u = vs(i) * 0.5 * parms.CB.s / parms.CB.h;

        if i == 1, X0 = parms.CB.Xmax;
        else,      X0 = Xmax(i-1,:);
        end

        % need to use fmincon because of nonlinear constraints violated with fminsearch
        [Xmax(i,:), cost(i)] = fmincon(@(X) cfxc.find_steadystate(X, parms), X0, [],[],[],[],LB(1:length(X0)),[],nlincon_func, opt);
        
        if cost(i) > 1e-3

            disp(['Could not solve concentric steady-state (',num2str(i),'/', num2str(length(vs)), '), simulating ...'])
            [~,x] = ode113(@cfxc.sim_muscle, [0 1], [1 parms.CB.Xmax], parms.set.odeopt, parms); 
            Xmax(i,:) = x(end,2:end);

        end
    end

    kCB_con = flip(Xmax(:,1)) / parms.CB.Xmax(1);
    FCB_con = flip(Xmax(:,2)) / parms.CB.Xmax(2);

    % eccentric
    id = v>0;
    vs = v(id);

    Xmax = nan(length(vs), length(X0));
    cost = nan(length(vs), 1);
    for i = 1:length(vs)
        parms.exp.u = vs(i) * 0.5 * parms.CB.s / parms.CB.h;

        if i == 1, X0 = parms.CB.Xmax;
        else,      X0 = Xmax(i-1,:);
        end

        % need to use fmincon because of nonlinear constraints violated with fminsearch
        [Xmax(i,:), cost(i)] = fmincon(@(X) cfxc.find_steadystate(X, parms), X0, [],[],[],[],LB(1:length(X0)),[],nlincon_func, opt);
        
        % if cost too large, we need to simulate
        if cost(i) > 1e-3
            
            disp(['Could not solve eccentric steady-state (',num2str(i),'/', num2str(length(vs)), '), simulating ...'])
           [~,x] = ode113(@cfxc.sim_muscle, [0 1], [1 parms.CB.Xmax], parms.set.odeopt, parms); 
           Xmax(i,:) = x(end,2:end);
            
        end
    end

    kCB_ecc = Xmax(:,1) / parms.CB.Xmax(1);
    FCB_ecc = Xmax(:,2) / parms.CB.Xmax(2);

    % combine concentric and eccentric
    kCB = [kCB_con; kCB_ecc];
    FCB = [FCB_con; FCB_ecc];
end

function[] = compare_fv(fv, parms)

        % simulate isometric
        parms.exp.u = 0;
        parms.exp.A = 1;
        parms.exp.a = 1;
        parms.type = 'crossbridge';
        
        [t,x] = ode113(@cfxc.dXfunc_shell, [0 .1], [0 0 0], [], parms);
        n = cfxc.n_func(parms.CB.Xmax(1), parms.CB.Xmax(2), parms.CB.Xmax(3), parms.CB);

        subplot(221);
        plot(parms.CB.xi, parms.CB.f_func(parms)); hold on
        plot(parms.CB.xi, parms.CB.g_func(parms),'--');
        xlabel('Strain (h)'); ylabel('Rate (s^{-1})')
        legend('Attachment','Detachment','location','best'); 
        legend boxoff

        subplot(222);
        plot(parms.CB.xi, n); hold on
        xlabel('Strain (h)'); ylabel('n (a.u.)');

        subplot(223);
        plot(t,x); hold on
        plot(t(end), parms.CB.Xmax, 'o','color',[.5 .5 .5], 'markerfacecolor', [.5 .5 .5])
        xlabel('Time (s)'); ylabel('Q (a.u.)')
        legend('Q_0','Q_1','Q_2','location','best'); 
        legend boxoff
        
        subplot(224)
        plot(fv.vHill, fv.FHill); hold on
        plot(fv.vHill, fv.FCB);
        xlim([-1 1/2]*parms.ce.vmaxrel)
        xlabel('Velocity (l_{opt}/s)'); 
        ylabel('Force (F_{max})')
        legend('Hill','CB','CB (2)','DM','location','best'); 
        legend boxoff
        box off

        titles = {'Rate functions', 'Isometric distribution', 'Isometric contraction','Force-velocity'};
        for i = 1:4
        subplot(2,2,i); 
        box off;
        title(titles{i})
        end
        title('Force-velocity')
        set(gcf,'units','normalized','position', [0 .3 .6 .4])


end

function[parms] = gen_funcs(parms)
    % crossbridge rate functions
    parms.CB.f_func = @(parms) parms.CB.f(1) .* (parms.CB.xi>0 & parms.CB.xi<=1) .* parms.CB.xi;
    
    parms.CB.g_func = @(parms) parms.CB.g(1) .* (parms.CB.xi >= 0) .* parms.CB.xi + ...
                           parms.CB.g(2) .* (parms.CB.xi < 0) + ...
                           parms.CB.g(3) .* (parms.CB.xi >= 1) .* (parms.CB.xi - 1);
                       
                       
    parms.CB.scale_rates = @(u, rate, mu) rate * (mu + u*(1-mu)); 
        
    % CE force-length
    parms.func.fce = @(L, parms) (L >= ((1-parms.ce.w) * parms.ce.lceopt) & L <= ((1+parms.ce.w) * parms.ce.lceopt))...
                     .* (parms.ce.c .* (L/parms.ce.lceopt).^2  - 2*parms.ce.c .* (L/parms.ce.lceopt) + parms.ce.c + 1);
                 
    parms.func.fv = @(a,F,Fisom,parms) -((Fisom + parms.ce.Arel.*a.^-.3) .* (parms.ce.vmaxrel*parms.ce.Arel) ./ ((F./a) + parms.ce.Arel.*a.^-.3) - (parms.ce.vmaxrel*parms.ce.Arel)) .* ...
           ((-((Fisom + parms.ce.Arel.*a.^-.3) .* (parms.ce.vmaxrel*parms.ce.Arel) ./ ((F./a) + parms.ce.Arel.*a.^-.3) - (parms.ce.vmaxrel*parms.ce.Arel))) < 0) + ...
            -(((parms.ce.vmaxrel*parms.ce.Arel) .* (Fisom + (-Fisom .* parms.ce.Fasymp))^2 / ((Fisom + parms.ce.Arel.*a.^-.3) .* parms.ce.Slopefactor)) ./ ((F./a) + (-Fisom .* parms.ce.Fasymp)) - ...
                           ((parms.ce.vmaxrel*parms.ce.Arel) .* (Fisom + (-Fisom .* parms.ce.Fasymp))^2 / ((Fisom + parms.ce.Arel.*a.^-.3) .* parms.ce.Slopefactor)) / (Fisom + (-Fisom .* parms.ce.Fasymp))) .* ...
                       ((-((Fisom + parms.ce.Arel.*a.^-.3) .* (parms.ce.vmaxrel*parms.ce.Arel) ./ ((F./a) + parms.ce.Arel.*a.^-.3) - (parms.ce.vmaxrel*parms.ce.Arel))) > 0);

    % SE force-length
    parms.func.fse = @(dlse_rel, parms) (dlse_rel<parms.see.sexm & dlse_rel > 0) .* (parms.see.sefm * 1/(exp(parms.see.sesh)-1) * (exp(parms.see.sesh/parms.see.sexm * dlse_rel)-1)) ...
                                      +(dlse_rel>=parms.see.sexm) .* (parms.see.sefm + parms.see.sesl * (dlse_rel - parms.see.sexm));

    parms.func.kse = @(dlse_rel, parms) (dlse_rel<parms.see.sexm) .* (parms.see.sesh/parms.see.sexm) .* parms.see.sefm .* 1/(exp(parms.see.sesh)-1) .* (exp(parms.see.sesh/parms.see.sexm * dlse_rel)) ...
                                      +(dlse_rel>=parms.see.sexm) .* parms.see.sesl;

    parms.func.lse = @(frel, parms) (frel < parms.see.sefm) .* (parms.see.sexm/parms.see.sesh * log(1 + frel*(exp(parms.see.sesh)-1)/parms.see.sefm)) ... 
                                   +(frel >=parms.see.sefm) .* (parms.see.sexm + (frel-parms.see.sefm)/parms.see.sesl);

    % parallel element
    parms.func.fpe = @(L,parms) ((parms.ce.Fmax / ((parms.pe.lpemax - parms.pe.lpe0) * parms.ce.lceopt).^2) * (L - (parms.pe.lpe0 * parms.ce.lceopt)).^2) .* (L > (parms.pe.lpe0 * parms.ce.lceopt));
    parms.func.kpe = @(L,parms) (2 * (parms.ce.Fmax / ((parms.pe.lpemax - parms.pe.lpe0) * parms.ce.lceopt).^2) * (L - (parms.pe.lpe0 * parms.ce.lceopt))) .* (L > (parms.pe.lpe0 * parms.ce.lceopt));

    % MTC
    parms.func.lmtc = @(phi, parms) parms.mtc.lmtc0 + parms.mtc.r * phi/180*pi;
    
    % CE dynamics
    parms.func.act = @(u,Ca,parms)(u-Ca)/(parms.ce.tau(1)*(u>Ca)+(parms.ce.tau(2)*(u<=Ca)));
    parms.func.fac =  @(Ca,r,parms)(Ca-r)/(parms.ce.tauR(1)*(Ca>r)+(parms.ce.tauR(2)*(Ca<=r)));
    
    % default input functions
    parms.exp.u_func = @(t, parms) parms.exp.A * (t < parms.exp.tstop);
end
  
function[parms] = gen_parms(parms)

    % crossbridge
    parms.CB.analytical = 1;
    parms.CB.f = 150; % initial value, will be fit
    parms.CB.g = [150 1000 150];  % initial value, will be fit
    parms.CB.mu = 1/3;
    parms.CB.xi = linspace(-10,10,1000);
    
    % only for original Huxley-type model
    parms.CB.nbins = 20000;
    parms.CB.xi0 = linspace(-50,50,parms.CB.nbins);
    
    % scaling
    parms.CB.h = 12*10^-9; % [m], crossbridge reach
    parms.CB.s = 2.64 * 10^-6;  % [m], sarcomere length
    parms.CB.Xmax = [1/2 1/4 1/6];
    
    % CE force-length
    parms.ce.w = .56;
    parms.ce.c = -1/parms.ce.w^2;
    parms.ce.amin = 1e-3; % neccesary for Hill-type model
    
    % CE force-velocity             
    parms.ce.Slopefactor = 2; % [] slope at v = 0
    parms.ce.Fasymp = 1.5; % [Fmax] force asymptote at v > 0

    % PE
    parms.pe.lpe0 = 1.0; % previously 1.2
    parms.pe.lpemax = 1.6;

    % pennation
    parms.ce.pennation = 0;
end
  
function[T,X] = stretch_protocol(as, vs, ts, parms, type)
    
    T = 0;
    X = parms.exp.x0;
    
    for i = 1:length(vs)
    
        if strcmp(type,'muscle')
            parms.exp.u = vs(i);
        elseif strcmp(type,'muscle-tendon')
            parms.exp.vmtc = vs(i);
        end
        
        parms.exp.A = as(i);
        tmax = ts(i);
        
        [t,x] = ode113(@cfxc.sim_muscle, [0 tmax], X(end,:), parms.set.odeopt, parms);

        T = [T; t+T(end)];
        X = [X; x];
    end
end

function[] = state_plot(t,x,titles, ls, color)

    S = size(x,2);

    if S < 6
        N = 1;
        M = 6;
    elseif S == 6
        N = 2;
        M = 3;
    elseif S > 6 && S < 9
        N = 2;
        M = 4;
    elseif S == 9
        N = 3;
        M = 3;
    elseif S > 9 && S < 12
        N = 3;
        M = 4;
    elseif S == 12
        N = 2;
        M = 6;
    elseif S > 12 && S < 15
        N = 2;
        M = 7;
    elseif S == 15
        N = 3;
        M = 5;
     elseif S == 16
        N = 4;
        M = 4;
     elseif S > 16
        N = 3;
        M = 6;
     end
    

    if nargin < 4
        ls = '-';
    end

    if nargin < 5
        color = 'blue';
    end

    for i = 1:S
        subplot(N,M,i);
        plot(t, x(:,i),'linestyle',ls,'color',color); hold on; box off

        if nargin > 2
            title(titles{i})
        end

    end

end
end
end