function [y_, y1st, y2nd, y3rd]=simult_FGRU(y0,dr,ex_,iorder,y1st_start,u_1_start)
% Simulates the model using a perturbation approach, given the path for the exogenous variables and the
% decision rules. This file builds on the simult_.m of Dynare 4.4 and has been adapted from the Dynare version to
% reflect the pruning and simulation approach of FGRU. It is hard-wired to
% the model under consideration and cannot be used for other models without
% adjustment. The specific parts begin at line 108
%
% INPUTS
%    y0         [double]   n*1 vector, initial value (n is the number of declared endogenous variables plus the number 
%                          of auxilliary variables for lags and leads)
%    dr         [struct]   matlab's structure where the reduced form solution of the model is stored.
%    ex_        [double]   T*q matrix of innovations.
%    iorder     [integer]  order of the taylor approximation.
%    y1st_start [double]   n*1 vector, initial value for linear term in
%                          pruned state space
%    u_1_start  [double]   n_exo*1 vector, initial value for linear shock term 
% 
% OUTPUTS
%    y_         [double]   n*(T+1) time series for the endogenous variables.
%
%    y1st       [double]   n*(T+1) time series for the linear term of the
%                          pruned state space
%    y2nd       [double]   n*(T+1) time series for the quadratic term of the
%                          pruned state space
%    y3rd       [double]   n*(T+1) time series for the cubic term of the
%                          pruned state space
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2014 Dynare Team and Johannes Pfeifer
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare and this software are distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% A copy of the GNU General Public License can be found at <http://www.gnu.org/licenses/>.

global M_ options_

iter = size(ex_,1);
endo_nbr = M_.endo_nbr;
exo_nbr = M_.exo_nbr;

y_ = zeros(size(y0,1),iter+M_.maximum_lag);
y_(:,1) = y0;
if nargout>1 
    y1st=zeros(size(y0,1),iter+M_.maximum_lag);
    y2nd=zeros(size(y0,1),iter+M_.maximum_lag); 
    y3rd=zeros(size(y0,1),iter+M_.maximum_lag);
end

% stoch_simul sets k_order_solver=1 if order=3, but does so only locally, so we
% have to do it here also
if options_.order == 3
    options_.k_order_solver = 1;
end

if ~options_.k_order_solver || (options_.k_order_solver && options_.pruning) %if k_order_pert is not used or if we do not use Dynare++ with k_order_pert
    if iorder==1
        y_(:,1) = y_(:,1)-dr.ys;
    end
end

if options_.k_order_solver && ~options_.pruning % Call dynare++ routines.
    ex_ = [zeros(M_.maximum_lag,M_.exo_nbr); ex_];
    switch options_.order
      case 1
        [err, y_] = dynare_simul_(1,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,exo_nbr, ...
                                  y_(dr.order_var,1),ex_',M_.Sigma_e,options_.DynareRandomStreams.seed,dr.ys(dr.order_var),...
                                  zeros(endo_nbr,1),dr.g_1);
      case 2
        [err, y_] = dynare_simul_(2,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,exo_nbr, ...
                                  y_(dr.order_var,1),ex_',M_.Sigma_e,options_.DynareRandomStreams.seed,dr.ys(dr.order_var),dr.g_0, ...
                                  dr.g_1,dr.g_2);
      case 3
        [err, y_] = dynare_simul_(3,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,exo_nbr, ...
                                  y_(dr.order_var,1),ex_',M_.Sigma_e,options_.DynareRandomStreams.seed,dr.ys(dr.order_var),dr.g_0, ...
                                  dr.g_1,dr.g_2,dr.g_3);
      otherwise
        error(['order = ' int2str(order) ' isn''t supported'])
    end
    mexErrCheck('dynare_simul_', err);
    y_(dr.order_var,:) = y_;
else
    if options_.block
        if M_.maximum_lag > 0
            k2 = dr.state_var;
        else
            k2 = [];
        end;
        order_var = 1:endo_nbr;
        dr.order_var = order_var;
    else
        k2 = dr.kstate(find(dr.kstate(:,2) <= M_.maximum_lag+1),[1 2]);
        k2 = k2(:,1)+(M_.maximum_lag+1-k2(:,2))*endo_nbr;
        order_var = dr.order_var;
    end;
    
    switch iorder
      case 3 % do FGRU pruning
        % only with pruning
        ghx = dr.ghx;
        ghu = dr.ghu;
        ghxx = dr.ghxx;
        ghxu = dr.ghxu;
        ghuu = dr.ghuu;
        ghs2 = dr.ghs2;
        ghxxx = dr.ghxxx;
        ghxxu = dr.ghxxu;
        ghxuu = dr.ghxuu;
        ghuuu = dr.ghuuu;
        ghxss = dr.ghxss;
        ghuss = dr.ghuss;
        threads = options_.threads.kronecker.A_times_B_kronecker_C;
        nspred = M_.nspred;
        ipred = M_.nstatic+(1:nspred);
        if nargin<5
            yhat1 = zeros(length(k2),1);
            yhat3 = y0(order_var(k2))-dr.ys(order_var(k2));
        elseif nargin>=5
            yhat1=y1st_start(order_var(k2));
            yhat3 = y0(order_var(k2))-dr.ys(order_var(k2));
        end
        exogenous_states={'sigma_r';          
            'sigma_tb'
            'eps_r'
            'eps_tb'
            'X'};
        ex_state_pos_dr=[]; %position in variable vector in decision rule order (LHS of equations below)
        for ii=1:size(exogenous_states,1)
            ex_state_pos_dr=[ex_state_pos_dr; strmatch(deblank(exogenous_states(ii,:)),M_.endo_names(order_var,:),'exact')];
        end
        for i=2:iter+M_.maximum_lag
            u = ex_(i-1,:)';
            %% setting different shock for first order term; reflects stale term from previous run in FGRU 
            if i==2
                u_1=u_1_start;
            else
                u_1=u;
            end
            %% construct terms of order 2 from second order part, based on linear component yhat1
            [gyy, err] = A_times_B_kronecker_C(ghxx,yhat1,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [guu, err] = A_times_B_kronecker_C(ghuu,u_1,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [gyu, err] = A_times_B_kronecker_C(ghxu,yhat1,u_1,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            %% construct terms of order 3, all based on first order component yhat1              
            y2a = kron(yhat1,yhat1);
            [gyyy, err] = A_times_B_kronecker_C(ghxxx,y2a,yhat1,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            u2a = kron(u_1,u_1);
            [guuu, err] = A_times_B_kronecker_C(ghuuu,u2a,u_1,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            yu = kron(yhat1,u_1);
            [gyyu, err] = A_times_B_kronecker_C(ghxxu,yhat1,yu,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [gyuu, err] = A_times_B_kronecker_C(ghxuu,yu,u_1,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            
            %% get unpruned solution for exogenous states 
            [gyy_unpruned, err] = A_times_B_kronecker_C(ghxx,yhat3,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [guu_unpruned, err] = A_times_B_kronecker_C(ghuu,u,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [gyu_unpruned, err] = A_times_B_kronecker_C(ghxu,yhat3,u,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            %% construct terms of order 3
            y2a_unpruned = kron(yhat3,yhat3);
            [gyyy_unpruned, err] = A_times_B_kronecker_C(ghxxx,y2a_unpruned,yhat3,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            u2a_unpruned = kron(u,u);
            [guuu_unpruned, err] = A_times_B_kronecker_C(ghuuu,u2a_unpruned,u,threads);
            yu_unpruned = kron(yhat3,u);
            [gyyu_unpruned, err] = A_times_B_kronecker_C(ghxxu,yhat3,yu_unpruned,threads);
            mexErrCheck('A_times_B_kronecker_C', err);
            [gyuu_unpruned, err] = A_times_B_kronecker_C(ghxuu,yu_unpruned,u,threads);
            mexErrCheck('A_times_B_kronecker_C', err);

            y_exo = ghx*yhat3 +ghu*u ...
                        + 1/2*(gyy_unpruned + guu_unpruned + 2*gyu_unpruned + ghs2) ...
                        + 1/6*(gyyy_unpruned + guuu_unpruned + 3*(gyyu_unpruned + gyuu_unpruned +  ghxss*yhat3 + ghuss*u)); %note: s is treated as variable, thus xss and uss are third order
            
            %% add all terms of order 3, linear component based on third order yhat3
            yhat3 = ghx*yhat3 +ghu*u ...
                    + 1/2*(gyy + guu + 2*gyu + ghs2) ...
                    + 1/6*(gyyy + guuu + 3*(gyyu + gyuu +  ghxss*yhat1 + ghuss*u_1)); %note: s is treated as variable, thus xss and uss are third order
            yhat1 = ghx*yhat1 + ghu*u_1;

            yhat3(ex_state_pos_dr,:)=y_exo(ex_state_pos_dr,:);
            yhat1(ex_state_pos_dr,:)=y_exo(ex_state_pos_dr,:);
            y_(order_var,i) = dr.ys(order_var)+ yhat3; %combine terms again
            y1st(order_var,i)= yhat1;
            y3rd(order_var,i)= yhat3;
            yhat1 = yhat1(ipred);
            yhat3 = yhat3(ipred);
        end
        y_(~ismember(1:M_.endo_nbr,dr.inv_order_var(k2)),1)=dr.ys(~ismember(1:M_.endo_nbr,dr.inv_order_var(k2))); % set first period that is not simulated to steady state
    end
end
