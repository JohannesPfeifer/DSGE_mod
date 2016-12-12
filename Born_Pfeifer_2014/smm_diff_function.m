function [fval, simulated_moments]=smm_diff_function(xopt,target,replications,shocks)
% function [fval, simulated_moments]=smm_diff_function(xopt,target,replications,shocks)
% Computes the quadratic deviation of the simulated moments from the target
% moments in the data
% Inputs:
%   xopt            [npar by 1]     vector of parameters in current
%                                   optimization step
%   target          [npar by 1]     vector of target data moments
%   replications    [scalar]        number of simulation replications
%   shocks          [nperiods by nshocks by replications] matrix of shocks
%                                   for simulation
% Outputs:
%   fval            [scalar]        value of the objective function (quadratic distance)
%   simulated_moments [npar by 1]   simulated moments in current step
% 
% The file is used to replicate the results of
% Benjamin Born and Johannes Pfeifer (2014): "Risk Matters: A comment", American Economic Review
% 
% Copyright (C) 2013-14 Benjamin Born and Johannes Pfeifer
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.

global oo_ M_ options_ % get Dynare structures; used to pass them on to resol.m
tic

%% simulate data from monthly model
npts=size(shocks,1);
out_withshock=NaN(M_.endo_nbr,npts,replications); %initialize matrix

%% set parameter for use in Dynare
set_param_value('sigma_x',xopt(1)); % TFP volatility
set_param_value('phipar',xopt(2)); % capital adjustment costs
set_param_value('D_bar',xopt(3));  % steady state debt
set_param_value('Phi',log(xopt(4))); %debt elasticity

if xopt(2)<0 %make sure adjustment cost is positive; using log-transformation as for Phi would yield numerical overflow
    fval=10e6+sum([xopt(1),xopt(2),xopt(3),log(xopt(4))].^2); %penalty function
    toc
    return
end

[oo_.dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_); %run model solution in Dynare

if info %solution was not successful
    fval=10e6+sum([xopt(1),xopt(2),xopt(3),log(xopt(4)) ].^2); %return with penalty 
else
    %% compute EMAS to get NX share
    burnin=4000;
    shock_mat=zeros(burnin,M_.exo_nbr);
    out_noshock = simult_FGRU(oo_.dr.ys,oo_.dr,shock_mat,options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1));
    ergodicmean_no_shocks=out_noshock(:,end);
    simulated_moments(1,4)=ergodicmean_no_shocks(strmatch('NX_Y',M_.endo_names,'exact'),1)*100; %net export share at EMAS

    %%simulate time series
    for ii=1:replications
        if ii==1
            [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shocks(1:end-1,:,ii),options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1)); %last shock is only needed as initial condition for next run
        else
            u_1_start=shocks(end,:,ii-1); %stale first order shock term from previous run
            [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shocks(1:end-1,:,ii),options_.order,y1st(:,end),u_1_start');        
        end
        out_withshock(:,:,ii) =y_; %use all including leading steady state
    end

    %% perform aggregation variables

    n_quarters=floor(npts/3);

    %initialize data strucuture for aggregated variables
    Y_quarterly_mean_aggregation=NaN(replications,n_quarters);
    C_quarterly_mean_aggregation=NaN(replications,n_quarters);
    I_quarterly_mean_aggregation=NaN(replications,n_quarters);
    NX_quarterly_mean_aggregation=NaN(replications,n_quarters);

    for ii=1:n_quarters
        %aggregation by mean over quarters
            Y_quarterly_mean_aggregation(:,ii)=squeeze(mean(out_withshock(strmatch('Y',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
            C_quarterly_mean_aggregation(:,ii)=squeeze(mean(out_withshock(strmatch('C',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
            I_quarterly_mean_aggregation(:,ii)=squeeze(mean(out_withshock(strmatch('I',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
            NX_quarterly_mean_aggregation(:,ii)=squeeze(mean(out_withshock(strmatch('NX',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));    
    end
    % get cyclical properties
    [~, Y_cyc_mean_aggregation]=hpfilter(Y_quarterly_mean_aggregation',1600);
    [~, C_cyc_mean_aggregation]=hpfilter(C_quarterly_mean_aggregation',1600);
    [~, I_cyc_mean_aggregation]=hpfilter(I_quarterly_mean_aggregation',1600);

    % compute target statistics
    simulated_moments(1,1)=mean(std(Y_cyc_mean_aggregation))*100;
    simulated_moments(1,2)=mean(std(C_cyc_mean_aggregation))*100/simulated_moments(1,1);
    simulated_moments(1,3)=mean(std(I_cyc_mean_aggregation))*100/simulated_moments(1,1);

    % compute objective function
    fval=(simulated_moments-target)*(simulated_moments-target)'
end
toc