% This file generates Appendix - Figure 2. Debt/Output, Current Account, and Net Exports Dynamicsof:
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

for country_iter=1:length(country_names)

    figure
    for recal_vs_replic_iter=1:2
            if recal_vs_replic_iter==1    
                load([country_names{1,country_iter},'_replication'],'ergodicmean_no_shocks','D_IRF_quarterly','Y_IRF_quarterly','CA_IRF_quarterly','NX_IRF_quarterly','M_','n_quarters');
            else
                load([country_names{1,country_iter},'_recalibration'],'ergodicmean_no_shocks','D_IRF_quarterly','Y_IRF_quarterly','CA_IRF_quarterly','NX_IRF_quarterly','M_','n_quarters');
            end
            %delt to output ratio
            subplot(3,2,recal_vs_replic_iter)
            D_Y_IRF=D_IRF_quarterly./Y_IRF_quarterly*100;
            plot(1:n_quarters,D_Y_IRF,'b-','LineWidth',1.5)
            title('D/Y','FontSize',10)
            axis tight
            set(gca,'FontSize',10)
            %current account to output ratio
            subplot(3,2,2+recal_vs_replic_iter)
            CA_Y_IRF_quarterly=CA_IRF_quarterly./Y_IRF_quarterly*100;
            plot(1:n_quarters,CA_Y_IRF_quarterly,1:n_quarters,zeros(1,n_quarters),'r--','LineWidth',1.5)
            title('CA/Y','FontSize',10)
            set(gca,'FontSize',10)
            axis tight
            %Net export to output ratio
            subplot(3,2,4+recal_vs_replic_iter)
            NX_Y_quarterly_mean=ergodicmean_no_shocks(strmatch('NX',M_.endo_names,'exact'),:)/exp(ergodicmean_no_shocks(strmatch('Y',M_.endo_names,'exact'),:))*100;
            NX_Y_percent_IRF_quarterly=NX_IRF_quarterly./Y_IRF_quarterly*100-NX_Y_quarterly_mean;
            plot(1:n_quarters,NX_Y_percent_IRF_quarterly,1:n_quarters,zeros(1,n_quarters),'r--','LineWidth',1.5)
            title('NX/Y','FontSize',10)
            axis tight
            set(gca,'FontSize',10)
    end
    print('-depsc2',['Figure_2_Appendix_D_NX_CA_response_',country_names{1,country_iter}])
end

