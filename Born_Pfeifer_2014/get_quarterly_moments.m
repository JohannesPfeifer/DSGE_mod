function [moments, row_label, column_label, std_NX]=get_quarterly_moments(simulated_series,EMAS,M_,oo_)
% function [moments, row_label, tb_y_corr,std_NX]=get_quarterly_moments(simulated_series,reference_point,M_,oo_)
% Aggregates the monthly to quarterly data and reports the moments for the 
% quarterly data using both the correct and the FGRU aggregation
% Inputs:
%   simulated_series    [nvars by nperiods] matrix containing the simulated
%                                           monthly series
%   EMAS                [nvars by 1] vector containing the EMAS, required
%                                       for reading out the net export to output ratio
%   M_                  [structure] Dynare model structure
%   oo_                 [structure] Dynare results structure
% 
% Outputs:
%   moments             [nmoments by 2] matrix storing moments; first
%                                       column: FGRU aggregation by sum; 
%                                       second column: correct aggregation
%                                       by mean
%   row_label           [nmoments by 1] string containing the moment names
%   column_label        [2 by 1]        string containing the column names
%   std_NX              [replication by 3] matrix with standard deviations
%                                       of the FGRU net export measure (first column), 
%                                       the net export to output ratio (second column), 
%                                       and the output volatility (third
%                                       column) over the simulation
%                                       repetitions
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

%% get info about dimensions
replications=size(simulated_series,3);
n_quarters=floor(size(simulated_series,2)/3);

%% initialize data strucuture for aggregated variables
Y_level_quarterly_sum_aggregation=NaN(replications,n_quarters);
% Correct aggregation
Y_quarterly_sum_aggregation=NaN(replications,n_quarters);
C_quarterly_sum_aggregation=NaN(replications,n_quarters);
I_quarterly_sum_aggregation=NaN(replications,n_quarters);
NX_quarterly_sum_aggregation=NaN(replications,n_quarters);
NX_Y_quarterly_sum_aggregation=NaN(replications,n_quarters);

% FGRU aggregation
Y_quarterly_mean_aggregation=NaN(replications,n_quarters);
C_quarterly_mean_aggregation=NaN(replications,n_quarters);
I_quarterly_mean_aggregation=NaN(replications,n_quarters);
NX_quarterly_mean_aggregation=NaN(replications,n_quarters);

%% Perform aggregation of variables

for ii=1:n_quarters
    %FGRU aggregation by sum over months
        Y_quarterly_sum_aggregation(:,ii)=squeeze(sum(simulated_series(strmatch('Y',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
        Y_level_quarterly_sum_aggregation(:,ii)=squeeze(sum(exp(simulated_series(strmatch('Y',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:))));
        C_quarterly_sum_aggregation(:,ii)=squeeze(sum(simulated_series(strmatch('C',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
        I_quarterly_sum_aggregation(:,ii)=squeeze(sum(simulated_series(strmatch('I',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
        NX_quarterly_sum_aggregation(:,ii)=squeeze(sum(simulated_series(strmatch('NX',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
        NX_Y_quarterly_sum_aggregation(:,ii)=squeeze(simulated_series(strmatch('NX_Y_quarterly',M_.endo_names,'exact'),ii*3,:)); %take value third month for quarter as it was defined backward-looking     
    %Correct aggregation by mean over months
        Y_quarterly_mean_aggregation(:,ii)=squeeze(mean(simulated_series(strmatch('Y',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
        C_quarterly_mean_aggregation(:,ii)=squeeze(mean(simulated_series(strmatch('C',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
        I_quarterly_mean_aggregation(:,ii)=squeeze(mean(simulated_series(strmatch('I',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));
        NX_quarterly_mean_aggregation(:,ii)=squeeze(mean(simulated_series(strmatch('NX',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:)));    
end
%compute correct net export/trade balance to output ratio
TB_Y_quarterly=NX_quarterly_sum_aggregation./Y_level_quarterly_sum_aggregation; % correct as ratio corrects for summation

%Define Net exports as in FGRU (i.e. wrong) and use Correia et al. (1995) transformation
net_exp_wrong_aggregation_wrong_definition = Y_quarterly_sum_aggregation-3*oo_.dr.ys(strmatch('Y',M_.endo_names,'exact')) - (C_quarterly_sum_aggregation-3*oo_.dr.ys(strmatch('C',M_.endo_names,'exact'))) - (I_quarterly_sum_aggregation-3*oo_.dr.ys(strmatch('I',M_.endo_names,'exact')));
net_exp_wrong_aggregation_wrong_definition = (net_exp_wrong_aggregation_wrong_definition./(abs(mean(net_exp_wrong_aggregation_wrong_definition,2))*ones(1,n_quarters)) - 1)/100 ; %additional division by 100 as in FGRU

%Do the same for correct computation of net exports, i.e. using the correct
%steady state values, but wrong aggregation
net_exp_mean_aggregation_correct_definition = Y_quarterly_mean_aggregation-oo_.dr.ys(strmatch('Y',M_.endo_names,'exact')) - (C_quarterly_mean_aggregation-oo_.dr.ys(strmatch('C',M_.endo_names,'exact'))) - (I_quarterly_mean_aggregation-oo_.dr.ys(strmatch('I',M_.endo_names,'exact')));
net_exp_mean_aggregation_correct_definition = (net_exp_mean_aggregation_correct_definition./(abs(mean(net_exp_mean_aggregation_correct_definition,2))*ones(1,n_quarters)) - 1)/100 ;

%Correia et al. (1995) transformation on model-consistent NX
NX_quarterly_sum_aggregation=(NX_quarterly_sum_aggregation./(abs(mean(NX_quarterly_sum_aggregation,2))*ones(1,n_quarters))-1)/100;  %additional division by 100 as in FV et al.
NX_quarterly_mean_aggregation=(NX_quarterly_mean_aggregation./(abs(mean(NX_quarterly_mean_aggregation,2))*ones(1,n_quarters))-1)/100;

%Compute Net Export Share in Ergodic Mean without shocks
nx_y_EMAS=EMAS(strmatch('NX_Y',M_.endo_names,'exact'),end);

%% compute cyclical components
%FGRU Aggregation
[~, Y_cyc_FGRU]=hpfilter(Y_quarterly_sum_aggregation',1600);
[~, C_cyc_FGRU]=hpfilter(C_quarterly_sum_aggregation',1600);
[~, I_cyc_FGRU]=hpfilter(I_quarterly_sum_aggregation',1600);
%Correct NX definition, but with incorrect aggegation as in FGRU
[~, NX_cyc_correct_definition_wrong_aggregation]=hpfilter(NX_quarterly_sum_aggregation',1600);
%NX definition via logs with incorrect aggegation as in FGRU
[~, NX_cyc_FGRU]=hpfilter(net_exp_wrong_aggregation_wrong_definition',1600);
%Cyclicality of Net-export-to-output ratio
[~, NX_Y_cyc]=hpfilter(NX_Y_quarterly_sum_aggregation',1600);


%Correct Aggregation
[~, Y_cyc_mean_aggregation]=hpfilter(Y_quarterly_mean_aggregation',1600);
[~, C_cyc_mean_aggregation]=hpfilter(C_quarterly_mean_aggregation',1600);
[~, I_cyc_mean_aggregation]=hpfilter(I_quarterly_mean_aggregation',1600);
%Correct NX definition and correct aggregation
[~, NX_cyc_mean_aggregation]=hpfilter(NX_quarterly_mean_aggregation',1600);
%Incorrect NX definition via logs, but with correct aggregation
[~, NX_cyc_FGRU_definition_mean_aggregation]=hpfilter(net_exp_mean_aggregation_correct_definition',1600);
%Cyclicality of Net-export-to-output ratio
[~, NX_Y_cyc_mean_aggregation]=hpfilter(TB_Y_quarterly',1600);

row_label={'$\sigma_Y$'
'$\sigma_C/\sigma_Y$' 
'$\sigma_I/\sigma_Y$' 
'$\sigma_{NX}/\sigma_Y$' 
'$\rho_{NX,Y}$' 
'$\sigma_{NX,log}/\sigma_Y$' 
'$\rho_{NX,log},Y$' 
'$\tilde NX/\tilde Y$'
'$\sigma_{NX/Y}$'
'$\rho{NX/Y,Y}$'};

column_label={'Incorrect_Aggregation','Correct_Aggregation'};
%% write FGRU aggregation moments to first column
moments(1,1)=mean(std(Y_cyc_FGRU))*100;
moments(2,1)=mean(std(C_cyc_FGRU))*100/moments(1,1);
moments(3,1)=mean(std(I_cyc_FGRU))*100/moments(1,1);
moments(4,1)=mean(std(NX_cyc_correct_definition_wrong_aggregation))*100/moments(1,1); 
moments(5,1)=mean(diag(corr(NX_cyc_correct_definition_wrong_aggregation,Y_cyc_FGRU)));
moments(6,1)=mean(std(NX_cyc_FGRU))*100/moments(1,1); 
moments(7,1)=mean(diag(corr(NX_cyc_FGRU,Y_cyc_FGRU)));
moments(8,1)=nx_y_EMAS*100;
moments(9,1)=mean(std(NX_Y_cyc))*100;
moments(10,1)=mean(diag(corr(NX_Y_cyc,Y_cyc_FGRU)));

%% write correct aggregation moments to first column
moments(1,2)=mean(std(Y_cyc_mean_aggregation))*100;
moments(2,2)=mean(std(C_cyc_mean_aggregation))*100/moments(1,2);
moments(3,2)=mean(std(I_cyc_mean_aggregation))*100/moments(1,2);
moments(4,2)=mean(std(NX_cyc_mean_aggregation))*100/moments(1,2);
moments(5,2)=mean(diag(corr(NX_cyc_mean_aggregation,Y_cyc_mean_aggregation)));
moments(6,2)=mean(std(NX_cyc_FGRU_definition_mean_aggregation))*100/moments(1,2);  
moments(7,2)=mean(diag(corr(NX_cyc_FGRU_definition_mean_aggregation,Y_cyc_mean_aggregation)));
moments(8,2)=nx_y_EMAS*100;
moments(9,2)=mean(std(NX_Y_cyc_mean_aggregation))*100;
moments(10,2)=mean(diag(corr(NX_Y_cyc_mean_aggregation,Y_cyc_mean_aggregation)));

%% save behavior of std(NX) over draws
std_NX=std(NX_cyc_mean_aggregation)';
std_NX=[std_NX,std(NX_Y_cyc_mean_aggregation)'];
std_NX=[std_NX,std(Y_cyc_mean_aggregation)'];
