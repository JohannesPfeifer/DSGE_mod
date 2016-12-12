% This file generates Table 2—Parameters Obtained by Moment Matching of:
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

fid=fopen('Table2_Parameter_table_Argentina.txt','w');

%% write header
header_string={'\begin{table}[!htbp]'
'  \centering'
'  \caption{Parameters Obtained by Moment Matching}'
'    \begin{tabularx}{0.80\textwidth}{r *{4}{Y}}'
'  \toprule'
'    & $\Phi_D$ & $\bar D$ &$\phi$ & $\sigma_x$ \\'
'    \cmidrule{2-5}'};
for ii=1:size(header_string,1)
   fprintf(fid,'%-120s\n',header_string{ii,1});
end


%load FGRU values from replication file
load('Argentina_replication','sigma_x','Phi','D_bar','phipar');
Parameters(1:4,1)=[Phi,D_bar,phipar,exp(sigma_x)]';

%load values from recalibration file
load('Argentina_recalibration','sigma_x','Phi','D_bar','phipar');
Parameters(1:4,2)=[Phi,D_bar,phipar,exp(sigma_x)]';

fprintf(fid,'%-30s \t & \t %5.2e \t & \t %5.2f \t & \t %5.2f \t & \t %5.3f \\\\ \n','Recalibration',Parameters(:,2));
fprintf(fid,'%-30s \t & \t %5.2e \t & \t %5.2f \t & \t %5.2f \t & \t %5.3f \\\\ \n','FGRU',Parameters(:,1));
    
    
bottom_string={'    \bottomrule'
'    \end{tabularx}'
'    \begin{tablenotes}'
'    first row: parameters obtained by moment matching using the corrected model. Second row: parameters obtained by moment matching in FGRU'
'    \end{tablenotes}'
'  \label{tab:recalibration}%'
'\end{table}%'};
for ii=1:size(bottom_string,1)
   fprintf(fid,'%-120s\n',bottom_string{ii,1});
end

fclose(fid)

edit Table2_Parameter_table_Argentina.txt