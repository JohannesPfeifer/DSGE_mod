% This file generates Table 3—Variance Decomposition of:
% Benjamin Born and Johannes Pfeifer (2014): "Risk Matters: A comment", American Economic Review
% 
% Requirements: requires that the main mod-file has been run before in
% replications and calibration mode and that the results have been
% correctly saved. See the ReadMe-file. Also requires a path set to Dynare.
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

%% store the FGRU results for volatility shocks only for reference in the last column
Argentina=[0.16;0.77;3.09;4.16];
Ecuador=[0.17;0.61;3.20;1.75];
Venezuela=[0.01;0.03;0.18;0.34];
Brazil=[0.04;0.08;0.54;0.26];

FGRU_var_decomp_vola_only=[Argentina,Ecuador,Venezuela,Brazil];
%% set up the LaTeX-Table for printing
fid=fopen('Table3_var_decomp.txt','w');

header_string={'\begin{table}[!tbp]'
'  \centering'
'  \caption{Variance Decomposition and Net Export Cyclicality}'
'    \begin{tabularx}{\textwidth}{r *{9}{Y}}'
'    \toprule'
'                              &        Data  &  i) All Shocks &    ii) TFP Only  &  iii) w/o Vola  &   iv) Rate  Level &  v) w/o TFP  &  vi) Vola Only  & FGRU Vola Only\\'};
for ii=1:size(header_string,1)
   fprintf(fid,'%-120s\n',header_string{ii,1});
end

for kk=1:length(country_names)
load([country_names{1,kk},'_recalibration']); %use _replication to generate decomposition to get the FGRU numbers after running the main mod-file in replication mode
replications=200;    
%% get ergodic mean in the absence of shocks as baseline
burnin=4000;
npts=96;
shock_mat=zeros(burnin,M_.exo_nbr);
out_noshock = simult_FGRU(oo_.dr.ys,oo_.dr,shock_mat,options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1));
log_deviations_SS_noshock=out_noshock -oo_.dr.ys*ones(1,burnin+M_.maximum_lag);
ergodicmean_no_shocks=out_noshock(:,end);

%% use FGRU shocks
winsorizing_dummy=1; %enable winsorizing
[shock_mat]=generate_FGRU_shocks(replications,winsorizing_dummy);
shock_mat_full=shock_mat;

%% baseline all shocks
shock_mat=shock_mat_full;

out_withshock=NaN(M_.endo_nbr,npts,replications);
for ii=1:replications
    %original calibration
    % use FV et al. pruning and start at EMAS
    if ii==1
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1)); %last shock is only needed as initial condition for next run
    else
        u_1_start=shock_mat(end,:,ii-1); %stale first order shock term from previous run
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,y1st(:,end),u_1_start');        
    end
    out_withshock(:,:,ii) =y_; %use all including leading steady state
end
[moments, column_label]=get_quarterly_moments(out_withshock,out_noshock,M_,oo_);
Var_decomp(:,1)=moments(:,2);

%% TFP only
shock_mat=zeros(npts,M_.exo_nbr,replications);
shock_mat(:,1,:)=shock_mat_full(:,1,:);
out_withshock=NaN(M_.endo_nbr,npts,replications);
for ii=1:replications
    %original calibration
    % use FV et al. pruning and start at EMAS
    if ii==1
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1)); %last shock is only needed as initial condition for next run
    else
        u_1_start=shock_mat(end,:,ii-1); %stale first order shock term from previous run
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,y1st(:,end),u_1_start');        
    end
    out_withshock(:,:,ii) =y_; %use all including leading steady state
end

[moments]=get_quarterly_moments(out_withshock,out_noshock,M_,oo_);
Var_decomp(:,2)=moments(:,2);

%% Without volatility
shock_mat=zeros(npts,M_.exo_nbr,replications);
shock_mat(:,1:3,:)=shock_mat_full(:,1:3,:);
out_withshock=NaN(M_.endo_nbr,npts,replications);
for ii=1:replications
    %original calibration
    % use FV et al. pruning and start at EMAS
    if ii==1
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1)); %last shock is only needed as initial condition for next run
    else
        u_1_start=shock_mat(end,:,ii-1); %stale first order shock term from previous run
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,y1st(:,end),u_1_start');        
    end
    out_withshock(:,:,ii) =y_; %use all including leading steady state
end
[moments]=get_quarterly_moments(out_withshock,out_noshock,M_,oo_);
Var_decomp(:,3)=moments(:,2);

%% Only Interest Rate Without volatility

shock_mat=zeros(npts,M_.exo_nbr,replications);
shock_mat(:,2:3,:)=shock_mat_full(:,2:3,:);
out_withshock=NaN(M_.endo_nbr,npts,replications);
for ii=1:replications
    %original calibration
    % use FV et al. pruning and start at EMAS
    if ii==1
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1)); %last shock is only needed as initial condition for next run
    else
        u_1_start=shock_mat(end,:,ii-1); %stale first order shock term from previous run
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,y1st(:,end),u_1_start');        
    end
    out_withshock(:,:,ii) =y_; %use all including leading steady state
end
[moments]=get_quarterly_moments(out_withshock,out_noshock,M_,oo_);
Var_decomp(:,4)=moments(:,2);

%% All without TFP

shock_mat=zeros(npts,M_.exo_nbr,replications);
shock_mat(:,2:5,:)=shock_mat_full(:,2:5,:);
out_withshock=NaN(M_.endo_nbr,npts,replications);
for ii=1:replications
    %original calibration
    % use FV et al. pruning and start at EMAS
    if ii==1
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1)); %last shock is only needed as initial condition for next run
    else
        u_1_start=shock_mat(end,:,ii-1); %stale first order shock term from previous run
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,y1st(:,end),u_1_start');        
    end
    out_withshock(:,:,ii) =y_; %use all including leading steady state
end
[moments]=get_quarterly_moments(out_withshock,out_noshock,M_,oo_);
Var_decomp(:,5)=moments(:,2);


%% All without TFP

shock_mat=zeros(npts,M_.exo_nbr,replications);
shock_mat(:,4:5,:)=shock_mat_full(:,4:5,:);
out_withshock=NaN(M_.endo_nbr,npts,replications);
for ii=1:replications
    %original calibration
    % use FV et al. pruning and start at EMAS
    if ii==1
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1)); %last shock is only needed as initial condition for next run
    else
        u_1_start=shock_mat(end,:,ii-1); %stale first order shock term from previous run
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,y1st(:,end),u_1_start');        
    end
    out_withshock(:,:,ii) =y_; %use all including leading steady state
end
[moments, header]=get_quarterly_moments(out_withshock,out_noshock,M_,oo_);
Var_decomp(:,6)=moments(:,2);


%% write everything to LaTeX-File
fprintf('\n%30s \n',country_names{1,kk}) 
fprintf('%25s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n','','Data','All Shocks','TFP Only','Without Vola','Rate Only','Without TFP','Volatility Only','FGRU Vola Only')
fprintf('%25s \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f\n','$\sigma_y$',moments_emp(1,1),Var_decomp(1,1),Var_decomp(1,2),Var_decomp(1,3),Var_decomp(1,4),Var_decomp(1,5),Var_decomp(1,6),FGRU_var_decomp_vola_only(1,kk))
header_string={'$\sigma_c$','$\sigma_i$','$\sigma_nx$'};
for ii=2:3
fprintf('%25s \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f\n',header_string{ii-1},moments_emp(ii,1)*moments_emp(1,1),Var_decomp(ii,1)*Var_decomp(1,1),Var_decomp(ii,2)*Var_decomp(1,2),Var_decomp(ii,3)*Var_decomp(1,3),Var_decomp(ii,4)*Var_decomp(1,4),Var_decomp(ii,5)*Var_decomp(1,5),Var_decomp(ii,6)*Var_decomp(1,6),FGRU_var_decomp_vola_only(ii,kk))
end
save([country_names{1,kk},'_',end_save_string,'_var_decomp_results'],'Var_decomp');

fprintf(fid,'%150s \n','    \midrule');
fprintf(fid,'%150s \n',['& \multicolumn{8}{c}{',country_names{1,kk},'}\\']); 
fprintf(fid,'%150s \n','    \cmidrule{2-9}');
fprintf(fid,'%25s \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \\\\\n','$\sigma_y$',moments_emp(1,1),Var_decomp(1,1),Var_decomp(1,2),Var_decomp(1,3),Var_decomp(1,4),Var_decomp(1,5),Var_decomp(1,6),FGRU_var_decomp_vola_only(1,kk));
for ii=2:3
fprintf(fid,'%25s \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \t & \t %10.2f \\\\\n',header_string{ii-1},moments_emp(ii,1)*moments_emp(1,1),Var_decomp(ii,1)*Var_decomp(1,1),Var_decomp(ii,2)*Var_decomp(1,2),Var_decomp(ii,3)*Var_decomp(1,3),Var_decomp(ii,4)*Var_decomp(1,4),Var_decomp(ii,5)*Var_decomp(1,5),Var_decomp(ii,6)*Var_decomp(1,6),FGRU_var_decomp_vola_only(ii,kk));
end

end

%% Write footer of Table
bottom_string={'    \bottomrule'
'    \end{tabularx}'
'    \begin{figurenotes}'
'            first column: moments obtained from HP-filtered data (1993Q1 - 2004Q3); second column: 200 simulations of the recalibrated model; third column: TFP shocks only; fourth column: without volatility shocks to spread and T-Bill rate; fifth column: only level shocks to the spread and the T-Bill rate; sixth column: without TFP shocks; seventh column: only shocks to the volatility of spreads and the T-Bill rate.'
'    \end{figurenotes}'
'  \label{tab:var_decomp_recalibration}%'
'\end{table}%'};
for ii=1:size(bottom_string,1)
   fprintf(fid,'%-120s\n',bottom_string{ii,1});
end
fclose(fid);

edit Table3_var_decomp.txt