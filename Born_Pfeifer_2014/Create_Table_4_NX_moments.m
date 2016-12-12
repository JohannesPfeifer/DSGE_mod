% This file generates Table 4—Cyclicality and Volatility of Net Exports of:
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
% country_string={'Argentina','Ecuador','Venezuela','Brazil'};
country_string={'Argentina'};

for country_iter=1:length(country_string)
    
    eval(['fid=fopen(''Table4_NX_moments_',country_string{1,country_iter},'.txt'',''w'');'])

    %% write table header
    header_string={'\begin{table}[htbp]'
    '  \centering'
    '  \caption{Cyclicality and Volatility of Net Exports}'
    '    \begin{tabularx}{\textwidth}{r *{6}{Y}}'
    '    \toprule'
    '          & \multicolumn{3}{c}{Orig.\ Calib.} & Recalib. & Data\\'
    '    \cmidrule(r{.75em}){2-4} \cmidrule(l{.75em}){5-5} \cmidrule(l{.75em}){6-6}'
    '          & FGRU  & TA & TA+NX & TA+NX  \\'
    '  \cmidrule{2-6}'};
    for ii=1:size(header_string,1)
       fprintf(fid,'%-120s\n',header_string{ii,1});
    end


    recal=load([country_string{1,country_iter},'_recalibration'],'column_names','row_names','moments_short','moments_emp','FGRU_moments');
    replic=load([country_string{1,country_iter},'_replication'],'column_names','row_names','moments_short','moments_emp','FGRU_moments');

    % NX-output correlation Correia et al. measure
    fprintf(fid,'%-30s \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f\\\\ \n','$\rho_{NX,Y}$',...        
        replic.FGRU_moments(strmatch('$\rho_{NX,log},Y$',replic.row_names,'exact'),:),...
        replic.moments_short(strmatch('$\rho_{NX,log},Y$',replic.row_names,'exact'),strmatch('Correct_Aggregation',replic.column_names,'exact')),...
        replic.moments_short(strmatch('$\rho_{NX,Y}$',replic.row_names,'exact'),strmatch('Correct_Aggregation',replic.column_names,'exact')),...
        recal.moments_short(strmatch('$\rho_{NX,Y}$',replic.row_names,'exact'),strmatch('Correct_Aggregation',recal.column_names,'exact')),...
        replic.moments_emp(strmatch('$\rho_{NX,Y}$',replic.row_names,'exact'),1));
    % Relative NX output volatility, Correia et al. measure
    fprintf(fid,'%-30s \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f\\\\ \n','$\sigma_{NX}/\sigma_Y$',...        
        replic.FGRU_moments(strmatch('$\sigma_{NX,log}/\sigma_Y$',replic.row_names,'exact'),:),...
        replic.moments_short(strmatch('$\sigma_{NX,log}/\sigma_Y$',replic.row_names,'exact'),strmatch('Correct_Aggregation',replic.column_names,'exact')),...
        replic.moments_short(strmatch('$\sigma_{NX}/\sigma_Y$',replic.row_names,'exact'),strmatch('Correct_Aggregation',replic.column_names,'exact')),...
        recal.moments_short(strmatch('$\sigma_{NX}/\sigma_Y$',replic.row_names,'exact'),strmatch('Correct_Aggregation',recal.column_names,'exact')),...
        replic.moments_emp(strmatch('$\sigma_{NX}/\sigma_Y$',replic.row_names,'exact'),1));
    % NX-output ratio correlation with output
    fprintf(fid,'%-30s \t & \t - \t & \t - \t & \t - \t & \t %5.2f \t & \t %5.2f\\\\ \n','$\rho_{NX/Y,Y}$',...        
        recal.moments_short(strmatch('$\rho{NX/Y,Y}$',replic.row_names,'exact'),strmatch('Correct_Aggregation',recal.column_names,'exact')),...
        replic.moments_emp(strmatch('$\rho{NX/Y,Y}$',replic.row_names,'exact'),1));
    % Volatility of NX to output ratio
    fprintf(fid,'%-30s \t & \t - \t & \t - \t & \t - \t & \t %5.2f \t & \t %5.2f\\\\ \n','$\sigma_{NX/Y}$',...        
        recal.moments_short(strmatch('$\sigma_{NX/Y}$',replic.row_names,'exact'),strmatch('Correct_Aggregation',recal.column_names,'exact')),...
        replic.moments_emp(strmatch('$\sigma_{NX/Y}$',replic.row_names,'exact'),1));
    
    %% write table footer
    bottom_string={'\bottomrule'
    '    \end{tabularx}'
    '    \begin{tablenotes}'
    '    first column: moments reported in FGRU. Second column: moments correcting the  time aggregation (TA). Third column: moments correcting the time aggregation and net export computation (TA+NX). Fourth column: moments obtained from the recalibrated corrected model. Fifth column: moments obtained from HP-filtered data. Simulations are conducted with 200 repetitions of 96 periods  using the same pruning, simulation, filtering, and winsorizing scheme as FGRU. Upper panel: based on CNR-approximation to net exports. Bottom panel: based on net-export-to-output ratio.'
    '\end{tablenotes}'
    '  \label{tab:Moments_NX}'
    '\end{table}'};
    for ii=1:size(bottom_string,1)
       fprintf(fid,'%-120s\n',bottom_string{ii,1});
    end
    fclose(fid);
    eval(['edit Table4_NX_moments_',country_string{1,country_iter},'.txt'])
end