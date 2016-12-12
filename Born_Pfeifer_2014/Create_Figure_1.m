% This file generates Figure 1. Comparison of original IRFs vs. IRFs of the
% recalibrated model of:
% Benjamin Born and Johannes Pfeifer (2014): "Risk Matters: A comment", American Economic Review
% 
% Requirements: requires that the main mod-file has been run before in
% replications and calibration mode and that the results have been
% correctly saved. See the ReadMe-file.
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

clear all
% country_names={'Argentina','Ecuador','Venezuela','Brazil'};
country_names={'Argentina'};

plot_vars={'Y';'C';'I';'H';'eps_r';'D'}; %Dynare name of variables to plot
plot_vars_heading={'Output';'Consumption';'Investment';'Hours';'Interest Rate Spread';'Debt'}; %Figure heading corresponding to variables

for country_iter=1:length(country_names)
    load([country_names{1,country_iter},'_replication'],'IRFS_quarterly_FGRU','IRFS_quarterly','M_','n_quarters'); %load previously saved results
    figure('name','Aggregation per Mean')
    for ii=1:6
        subplot(3,2,ii)
        load([country_names{1,country_iter},'_replication'],'IRFS_quarterly_FGRU');
        IRF_plot=IRFS_quarterly_FGRU(strmatch(plot_vars(ii,:),M_.endo_names,'exact'),:);
        load([country_names{1,country_iter},'_recalibration'],'IRFS_quarterly');
        IRF_plot_recal=IRFS_quarterly(strmatch(plot_vars(ii,:),M_.endo_names,'exact'),:);
        if strmatch(plot_vars(ii,:),'Y','exact')
            Y_IRF_plot(country_iter,:)=IRF_plot;
            Y_IRF_plot_recal(country_iter,:)=IRF_plot_recal;
        end
        if max(abs(IRF_plot))  >1e-10
            plot(1:n_quarters,IRF_plot,'b-',1:n_quarters,IRF_plot_recal,'r--','LineWidth',1.5)
        else
            plot(1:n_quarters,zeros(1,n_quarters),'b-',1:n_quarters,zeros(1,n_quarters),'r--','LineWidth',1.5) 
        end
        title(plot_vars_heading(ii,:))
        if ii==5
            legend('FGRU IRF','Recalibration')
        end
        axis tight
        grid on
        hold on
        plot((1:n_quarters),zeros(n_quarters,1),'r')
    end

    print('-depsc2',['Figure_1_IRF_comparison_',country_names{1,country_iter}])
end
