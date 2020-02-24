% This file replicates Table 5 (Welfare comparison with an inefficient steady 
% state) of Born/Pfeifer (2018): "The New Keynesian Wage Phillips Curve: Calvo vs. Rotemberg". 
% It requires Dynare 4.5 to be in the path.

% Copyright (C) 2018 Johannes Pfeifer and Benjamin Born
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
% For a copy of the GNU General Public License,
% see <http://www.gnu.org/licenses/>.

if ~isdir('Latex/Inefficient')
   mkdir('Latex/Inefficient') 
end

Ramsey=0; %do not do Ramsey, planner objective at order=2 not implemented

%% EHL Calvo
if Ramsey==1
    dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=0 -Dstrict_targeting=1 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=1 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
    pause(1);
    results_EHL_Calvo_mat(:,1)=100*[values_technology; values_demand];
end
% strict targeting
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=0 -Dstrict_targeting=1 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Calvo_mat(:,2)=100*[values_technology; values_demand];
!copy Born_Pfeifer_2018_welfare_original_content.tex "Latex/Inefficient/EHL_Calvo_model_equations.tex"
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=0 -Dstrict_targeting=1 -Dprice_targeting=0 -Dwage_targeting=1 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Calvo_mat(:,3)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=0 -Dstrict_targeting=1 -Dprice_targeting=0 -Dwage_targeting=0 -Dcomposite_targeting=1 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Calvo_mat(:,4)=100*[values_technology; values_demand];
%% EHL Calvo
% flexible targeting
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=0 -Dstrict_targeting=0 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Calvo_mat(:,5)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=0 -Dstrict_targeting=0 -Dprice_targeting=0 -Dwage_targeting=1 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Calvo_mat(:,6)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=0 -Dstrict_targeting=0 -Dprice_targeting=0 -Dwage_targeting=0 -Dcomposite_targeting=1 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Calvo_mat(:,7)=100*[values_technology; values_demand];

results_EHL_Calvo_mat_display=results_EHL_Calvo_mat;

options_.noprint=0;
headers_string={' ';'Ramsey';'Price';'Wage';'Composite';'Price';'Wage';'Composite'};
labels_string={'sigma(pi_p)';'sigma(pi_w)';'sigma(tilde y)';'W unc.';'W cond.';'sigma(pi_p)';'sigma(pi_w)';'sigma(tilde y)';'W unc.';'W cond.'};
labels_string_tex={'\sigma(\pi_p)';'\sigma(\pi_w)';'\sigma(\tilde y)';'\lambda_{unc}';'\lambda_{cond}';'\sigma(\pi_p)';'\sigma(\pi_w)';'\sigma(\tilde y)';'\lambda_{unc}';'\lambda_{cond}'};
dyntable(options_,'EHLCalvo',headers_string,labels_string,results_EHL_Calvo_mat_display,size(labels_string,2)+2,5,4)      
dyn_latex_table_modified(options_,'EHLCalvo','Latex/Inefficient/EHL_Calvo_welfare',headers_string,labels_string_tex,results_EHL_Calvo_mat_display,size(labels_string,2)+2,8,6);

clear results_EHL_Calvo_mat
%% EHL Rotemberg
if Ramsey==1
    dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=0 -Dstrict_targeting=1 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=1 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
    pause(1);
    results_EHL_Rotemberg_mat(:,1)=100*[values_technology; values_demand];
end
%strict targeting
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=0 -Dstrict_targeting=1 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Rotemberg_mat(:,2)=100*[values_technology; values_demand];
!copy Born_Pfeifer_2018_welfare_original_content.tex "Latex/Inefficient/EHL_Rotemberg_model_equations.tex"
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=0 -Dstrict_targeting=1 -Dprice_targeting=0 -Dwage_targeting=1 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Rotemberg_mat(:,3)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=0 -Dstrict_targeting=1 -Dprice_targeting=0 -Dwage_targeting=0 -Dcomposite_targeting=1 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Rotemberg_mat(:,4)=100*[values_technology; values_demand];
%% EHL Rotemberg
% flexible targeting
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=0 -Dstrict_targeting=0 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Rotemberg_mat(:,5)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=0 -Dstrict_targeting=0 -Dprice_targeting=0 -Dwage_targeting=1 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Rotemberg_mat(:,6)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=0 -Dstrict_targeting=0 -Dprice_targeting=0 -Dwage_targeting=0 -Dcomposite_targeting=1 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_EHL_Rotemberg_mat(:,7)=100*[values_technology; values_demand];

results_EHL_Rotemberg_mat_display=results_EHL_Rotemberg_mat;

dyntable(options_,'EHLRotemberg',headers_string,labels_string,results_EHL_Rotemberg_mat_display,size(labels_string,2)+2,5,4)      
dyn_latex_table_modified(options_,'EHLRotemberg','Latex/Inefficient/EHL_Rotemberg_welfare',headers_string,labels_string_tex,results_EHL_Rotemberg_mat_display,size(labels_string,2)+2,8,6);

clear results_EHL_Rotemberg_mat

%% SGU Calvo
if Ramsey==1
    dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=1 -Dstrict_targeting=1 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=1 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
    pause(1);
    results_SGU_Calvo_mat(:,1)=100*[values_technology; values_demand];
end
% strict targeting
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=1 -Dstrict_targeting=1 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Calvo_mat(:,2)=100*[values_technology; values_demand];
!copy Born_Pfeifer_2018_welfare_original_content.tex "Latex/Inefficient/SGU_Calvo_model_equations.tex"
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=1 -Dstrict_targeting=1 -Dprice_targeting=0 -Dwage_targeting=1 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Calvo_mat(:,3)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=1 -Dstrict_targeting=1 -Dprice_targeting=0 -Dwage_targeting=0 -Dcomposite_targeting=1 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Calvo_mat(:,4)=100*[values_technology; values_demand];

%% SGU Calvo
% flexible targeting
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=1 -Dstrict_targeting=0 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Calvo_mat(:,5)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=1 -Dstrict_targeting=0 -Dprice_targeting=0 -Dwage_targeting=1 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Calvo_mat(:,6)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=1 -DSGU_framework=1 -Dstrict_targeting=0 -Dprice_targeting=0 -Dwage_targeting=0 -Dcomposite_targeting=1 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Calvo_mat(:,7)=100*[values_technology; values_demand];

results_SGU_Calvo_mat_display=results_SGU_Calvo_mat;

dyntable(options_,'SGUCalvo',headers_string,labels_string,results_SGU_Calvo_mat_display,size(labels_string,2)+2,5,4)      
dyn_latex_table_modified(options_,'SGUCalvo','Latex/Inefficient/SGU_Calvo_welfare',headers_string,labels_string_tex,results_SGU_Calvo_mat_display,size(labels_string,2)+2,8,6);

clear results_SGU_Calvo_mat

%% SGU Rotemberg
if Ramsey==1
    dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=1 -Dstrict_targeting=1 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=1 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
    pause(1);
    results_SGU_Rotemberg_mat(:,1)=100*[values_technology; values_demand];
end
% strict targeting

dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=1 -Dstrict_targeting=1 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Rotemberg_mat(:,2)=100*[values_technology; values_demand];
!copy Born_Pfeifer_2018_welfare_original_content.tex "Latex/Inefficient/SGU_Rotemberg_model_equations.tex"
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=1 -Dstrict_targeting=1 -Dprice_targeting=0 -Dwage_targeting=1 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Rotemberg_mat(:,3)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=1 -Dstrict_targeting=1 -Dprice_targeting=0 -Dwage_targeting=0 -Dcomposite_targeting=1 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Rotemberg_mat(:,4)=100*[values_technology; values_demand];

%% SGU Rotemberg
% flexible targeting
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=1 -Dstrict_targeting=0 -Dprice_targeting=1 -Dwage_targeting=0 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Rotemberg_mat(:,5)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=1 -Dstrict_targeting=0 -Dprice_targeting=0 -Dwage_targeting=1 -Dcomposite_targeting=0 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Rotemberg_mat(:,6)=100*[values_technology; values_demand];
dynare Born_Pfeifer_2018_welfare -DCalvo=0 -DSGU_framework=1 -Dstrict_targeting=0 -Dprice_targeting=0 -Dwage_targeting=0 -Dcomposite_targeting=1 -DRamsey_policy=0 -Dfixed_WPC_slope=1 -Defficient_steady_state=0
pause(1);
results_SGU_Rotemberg_mat(:,7)=100*[values_technology; values_demand];

results_SGU_Rotemberg_mat_display=results_SGU_Rotemberg_mat;

dyntable(options_,'SGURotemberg',headers_string,labels_string,results_SGU_Rotemberg_mat_display,size(labels_string,2)+2,5,4)      
dyn_latex_table_modified(options_,'SGURotemberg','Latex/Inefficient/SGU_Rotemberg_welfare',headers_string,labels_string_tex,results_SGU_Rotemberg_mat_display,size(labels_string,2)+2,8,6);

clear results_SGU_Rotemberg_mat

fidTeX = fopen('Latex/Inefficient/Welfare_summary.tex','w');

if Ramsey==0
    
    header_string={'\begin{sidewaystable}'
        '    \npdecimalsign{.}'
        '    \nprounddigits{3}'
        '    \centering'
        '    \caption{Welfare: Inefficient Steady State}'
        '    \label{tab:welfare_inefficient}'
        '    \begin{tabular}{cn{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}}'
        '    \toprule'
        '    & \multicolumn{6}{c}{EHL Calvo} & \multicolumn{6}{c}{EHL Rotemberg} \\'
        '    \cmidrule(rl){2-7} \cmidrule(rl){8-13}'
        '    &\multicolumn{3}{c}{Strict Targeting} & \multicolumn{3}{c}{Flexible Targeting} &\multicolumn{3}{c}{Strict Targeting} & \multicolumn{3}{c}{Flexible Targeting} \\'
        '    \cmidrule(rl){2-4}\cmidrule(rl){5-7} \cmidrule(rl){8-10}\cmidrule(rl){11-13}'
        '                  	          & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.}	 & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.} & 	     	         {Price}	 & 	      {Wage}	 & 	     {Comp.}	 & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.}\\'
        '    \cmidrule(rl){2-7}\cmidrule(rl){8-13}'
        '   & \multicolumn{12}{c}{Technology Shock} \\'
        '    \cmidrule(rl){2-13}'
        };
    
    for ii=1:size(header_string,1)
        fprintf(fidTeX,'%s\n',header_string{ii});
    end
    
    label_format_leftbound  = sprintf('$%%-%ds$',size(labels_string_tex,2)+1);
    value_format='%8.6f';
    for ii=4:size(results_EHL_Calvo_mat_display,1)/2
        fprintf(fidTeX,label_format_leftbound,deblank(labels_string_tex{ii,:}));
        fprintf(fidTeX,['\t & \t' value_format],[results_EHL_Calvo_mat_display(ii,2:end) results_EHL_Rotemberg_mat_display(ii,2:end)]);
        fprintf(fidTeX,' \\\\ \n');
    end
    
    fprintf(fidTeX,'    \\cmidrule(rl){2-13}')
    fprintf(fidTeX,'   & \\multicolumn{12}{c}{Demand Shock} \\\\')
    fprintf(fidTeX,'    \\cmidrule(rl){2-13}')
    
    
    for ii=size(results_EHL_Calvo_mat_display,1)/2+4:size(results_EHL_Calvo_mat_display,1)
        fprintf(fidTeX,label_format_leftbound,deblank(labels_string_tex{ii,:}));
        fprintf(fidTeX,['\t & \t' value_format],[results_EHL_Calvo_mat_display(ii,2:end) results_EHL_Rotemberg_mat_display(ii,2:end)]);
        fprintf(fidTeX,' \\\\ \n');
    end
    
    
    middle_string={'    \midrule'
        '    & \multicolumn{6}{c}{SGU Calvo} & \multicolumn{6}{c}{SGU Rotemberg} \\'
        '    \cmidrule(rl){2-7} \cmidrule(rl){8-13}'
        '   & \multicolumn{12}{c}{Technology Shock} \\'
        '    \cmidrule(rl){2-13}'
        };
    
    for ii=1:size(middle_string,1)
        fprintf(fidTeX,'%s\n',middle_string{ii});
    end
    
    label_format_leftbound  = sprintf('$%%-%ds$',size(labels_string_tex,2)+1);
    value_format='%8.6f';
    for ii=4:size(results_EHL_Calvo_mat_display,1)/2
        fprintf(fidTeX,label_format_leftbound,deblank(labels_string_tex{ii,:}));
        fprintf(fidTeX,['\t & \t' value_format],[results_SGU_Calvo_mat_display(ii,2:end) results_SGU_Rotemberg_mat_display(ii,2:end)]);
        fprintf(fidTeX,' \\\\ \n');
    end
    
    fprintf(fidTeX,'    \\cmidrule(rl){2-13}')
    fprintf(fidTeX,'   & \\multicolumn{12}{c}{Demand Shock} \\\\')
    fprintf(fidTeX,'    \\cmidrule(rl){2-13}')
    
    
    for ii=size(results_EHL_Calvo_mat_display,1)/2+4:size(results_EHL_Calvo_mat_display,1)
        fprintf(fidTeX,label_format_leftbound,deblank(labels_string_tex{ii,:}));
        fprintf(fidTeX,['\t & \t' value_format],[results_SGU_Calvo_mat_display(ii,2:end) results_SGU_Rotemberg_mat_display(ii,2:end)]);
        fprintf(fidTeX,' \\\\ \n');
    end
    
else
    header_string={'\begin{sidewaystable}'
        '    \npdecimalsign{.}'
        '    \nprounddigits{3}'
        '    \centering'
        '    \caption{Welfare: Inefficient Steady State}'
        '    \label{tab:welfare_inefficient}'
        '    \begin{tabular}{cn{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}n{1}{3}}'
        '    \toprule'
        '    & \multicolumn{7}{c}{EHL Calvo} & \multicolumn{7}{c}{EHL Rotemberg} \\'
        '    \cmidrule(rl){2-8} \cmidrule(rl){9-15}'
        '    & {Optimal} &\multicolumn{3}{c}{Strict Targeting} & \multicolumn{3}{c}{Flexible Targeting} & {Optimal} &\multicolumn{3}{c}{Strict Targeting} & \multicolumn{3}{c}{Flexible Targeting} \\'
        '    \cmidrule(rl){3-5}\cmidrule(rl){6-8} \cmidrule(rl){10-12}\cmidrule(rl){13-15}'
        '                  	         & 	     	 & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.}	 & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.} & 	     	         & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.}	 & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.}\\'
        '    \cmidrule(rl){2-8}\cmidrule(rl){9-15}'
        '   & \multicolumn{14}{c}{Technology Shock} \\'
        '    \cmidrule(rl){2-15}'
        };
    
    for ii=1:size(header_string,1)
        fprintf(fidTeX,'%s\n',header_string{ii});
    end
    
    label_format_leftbound  = sprintf('$%%-%ds$',size(labels_string_tex,2)+1);
    value_format='%8.6f';
    for ii=1:size(results_EHL_Calvo_mat_display,1)/2
        fprintf(fidTeX,label_format_leftbound,deblank(labels_string_tex{ii,:}));
        fprintf(fidTeX,['\t & \t' value_format],[results_EHL_Calvo_mat_display(ii,:) results_EHL_Rotemberg_mat_display(ii,:)]);
        fprintf(fidTeX,' \\\\ \n');
    end
    
    fprintf(fidTeX,'    \\cmidrule(rl){2-15}')
    fprintf(fidTeX,'   & \\multicolumn{14}{c}{Demand Shock} \\\\')
    fprintf(fidTeX,'    \\cmidrule(rl){2-15}')
    
    for ii=size(results_EHL_Calvo_mat_display,1)/2+1:size(results_EHL_Calvo_mat_display,1)
        fprintf(fidTeX,label_format_leftbound,deblank(labels_string_tex{ii,:}));
        fprintf(fidTeX,['\t & \t' value_format],[results_EHL_Calvo_mat_display(ii,:) results_EHL_Rotemberg_mat_display(ii,:)]);
        fprintf(fidTeX,' \\\\ \n');
    end
    
    
    middle_string={'    \midrule'
        '    & \multicolumn{7}{c}{SGU Calvo} & \multicolumn{7}{c}{SGU Rotemberg} \\'
        '    \cmidrule(rl){2-8} \cmidrule(rl){9-15}'
        '    & {Optimal} &\multicolumn{3}{c}{Strict Targeting} & \multicolumn{3}{c}{Flexible Targeting} & {Optimal} &\multicolumn{3}{c}{Strict Targeting} & \multicolumn{3}{c}{Flexible Targeting} \\'
        '    \cmidrule(rl){3-5}\cmidrule(rl){6-8} \cmidrule(rl){10-12}\cmidrule(rl){13-15}'
        '                  	         & 	     	 & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.}	 & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.} & 	     	         & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.}	 & 	     {Price}	 & 	      {Wage}	 & 	     {Comp.}\\'
        '    \cmidrule(rl){2-8}\cmidrule(rl){9-15}'
        '   & \multicolumn{14}{c}{Technology Shock} \\'
        '    \cmidrule(rl){2-15}'
        };
    
    for ii=1:size(middle_string,1)
        fprintf(fidTeX,'%s\n',middle_string{ii});
    end
    
    label_format_leftbound  = sprintf('$%%-%ds$',size(labels_string_tex,2)+1);
    value_format='%8.6f';
    for ii=1:size(results_EHL_Calvo_mat_display,1)/2
        fprintf(fidTeX,label_format_leftbound,deblank(labels_string_tex{ii,:}));
        fprintf(fidTeX,['\t & \t' value_format],[results_SGU_Calvo_mat_display(ii,:) results_SGU_Rotemberg_mat_display(ii,:)]);
        fprintf(fidTeX,' \\\\ \n');
    end
    
    fprintf(fidTeX,'    \\cmidrule(rl){2-15}')
    fprintf(fidTeX,'   & \\multicolumn{14}{c}{Demand Shock} \\\\')
    fprintf(fidTeX,'    \\cmidrule(rl){2-15}')
    
    for ii=size(results_EHL_Calvo_mat_display,1)/2+1:size(results_EHL_Calvo_mat_display,1)
        fprintf(fidTeX,label_format_leftbound,deblank(labels_string_tex{ii,:}));
        fprintf(fidTeX,['\t & \t' value_format],[results_SGU_Calvo_mat_display(ii,:) results_SGU_Rotemberg_mat_display(ii,:)]);
        fprintf(fidTeX,' \\\\ \n');
    end
    
    
end

bottom_string={'    \bottomrule'
    '    \end{tabular}'
    '\end{sidewaystable}'};

for ii=1:size(bottom_string,1)
    fprintf(fidTeX,'%s\n',bottom_string{ii});
end
fclose(fidTeX);