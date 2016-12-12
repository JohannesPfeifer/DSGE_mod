% This file generates Table 1—Targeted Moments of:
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



% country_string={'Argentina','Ecuador','Venezuela','Brazil'};
clear all
country_string={'Argentina'};

for country_iter=1:length(country_string)
    
    eval(['fid=fopen(''Table1_moments_',country_string{1,country_iter},'.txt'',''w'');'])
    %% write header string
    header_string={' \begin{table}[!htbp]'
    '   \centering'
    '   \caption{Targeted Moments}'
    ' \begin{tabularx}{\textwidth}{r *{5}{Y}}'
    ' \toprule'
    '                    & $\sigma_Y$ & $\sigma_C/\sigma_Y$ & $\sigma_I/\sigma_Y$& $\widetilde{NX}/\widetilde{Y}$ \\'
    '                    \cmidrule{2-5}'};
    for ii=1:size(header_string,1)
       fprintf(fid,'%-120s\n',header_string{ii,1});
    end


    recal=load([country_string{1,country_iter},'_recalibration'],'column_names','row_names','moments_short','moments_emp','FGRU_moments');
    replic=load([country_string{1,country_iter},'_replication'],'column_names','row_names','moments_short','moments_emp','FGRU_moments');

    fprintf(fid,'%-30s \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t \\\\ \n','FGRU',replic.FGRU_moments(strmatch('$\sigma_Y$',replic.row_names,'exact'),:),replic.FGRU_moments(strmatch('$\sigma_C/\sigma_Y$',replic.row_names,'exact'),:),replic.FGRU_moments(strmatch('$\sigma_I/\sigma_Y$',replic.row_names,'exact'),:),replic.FGRU_moments(strmatch('$\tilde NX/\tilde Y$',replic.row_names,'exact'),:));
    fprintf(fid,'%-30s \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t \\\\ \n','Corr. Aggreg.',replic.moments_short(strmatch('$\sigma_Y$',replic.row_names,'exact'),strmatch('Correct_Aggregation',replic.column_names,'exact')),replic.moments_short(strmatch('$\sigma_C/\sigma_Y$',replic.row_names,'exact'),strmatch('Correct_Aggregation',replic.column_names,'exact')),replic.moments_short(strmatch('$\sigma_I/\sigma_Y$',replic.row_names,'exact'),strmatch('Correct_Aggregation',replic.column_names,'exact')),replic.moments_short(strmatch('$\tilde NX/\tilde Y$',replic.row_names,'exact'),strmatch('Correct_Aggregation',replic.column_names,'exact')));
    fprintf(fid,'%-30s \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t \\\\ \n','Recalibration',recal.moments_short(strmatch('$\sigma_Y$',recal.row_names,'exact'),strmatch('Correct_Aggregation',recal.column_names,'exact')),recal.moments_short(strmatch('$\sigma_C/\sigma_Y$',recal.row_names,'exact'),strmatch('Correct_Aggregation',recal.column_names,'exact')),recal.moments_short(strmatch('$\sigma_I/\sigma_Y$',recal.row_names,'exact'),strmatch('Correct_Aggregation',recal.column_names,'exact')),recal.moments_short(strmatch('$\tilde NX/\tilde Y$',recal.row_names,'exact'),strmatch('Correct_Aggregation',recal.column_names,'exact')));
    fprintf(fid,'%-30s \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t & \t %5.2f \t \\\\ \n','Data',recal.moments_emp(strmatch('$\sigma_Y$',replic.row_names,'exact'),1),recal.moments_emp(strmatch('$\sigma_C/\sigma_Y$',replic.row_names,'exact'),:),recal.moments_emp(strmatch('$\sigma_I/\sigma_Y$',replic.row_names,'exact'),:),recal.moments_emp(strmatch('$\tilde NX/\tilde Y$',replic.row_names,'exact'),:));
    %% write bottom string
    bottom_string={'    \bottomrule'
    '    \end{tabularx}'
    '    \begin{tablenotes}'
    'First row: moments reported in FGRU. Second row: FGRU with corrected time aggregation. Third row: moments obtained from simulating the recalibrated corrected model 200 times for 96 periods using the same pruning, simulation, filtering, and winsorizing scheme as FGRU. Fourth row: Moments obtained from HP-filtered data (1993Q1 - 2004Q3).'
    '\end{tablenotes}'
    '  \label{tab:Moments_recalibration}'
    '\end{table}'};
    for ii=1:size(bottom_string,1)
       fprintf(fid,'%-120s\n',bottom_string{ii,1});
    end

    fclose(fid);
    eval(['edit Table1_moments_',country_string{1,country_iter},'.txt'])
end
