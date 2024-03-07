clearvars all; clearvars -global;

%% Aguiar_Gopinath_2007
try
    close all;
    cd('Aguiar_Gopinath_2007');
    dynare Aguiar_Gopinath_2007
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Andreasen_2012
try
    close all; clearvars all; clearvars -global;
    cd('../Andreasen_2012');
    dynare Andreasen_2012_rare_disasters
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Ascari_Sbordone_2014
try
    close all; clearvars all; clearvars -global;
    cd('../Ascari_Sbordone_2014');
    dynare Ascari_Sbordone_2014
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Basu_Bundick_2017
try
    close all; clearvars all; clearvars -global;
    cd('../Basu_Bundick_2017');
    dynare Basu_Bundick_2017
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Born_Pfeifer_2014
try
    close all; clearvars all; clearvars -global;
    cd('../Born_Pfeifer_2014');
    dynare Born_Pfeifer_RM_Comment.mod
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Born_Pfeifer_2018
try
    close all; clearvars all; clearvars -global;
    cd('../Born_Pfeifer_2018/Monetary_Policy_IRFs');
    run_IRF_comparison;
    cd('../Welfare');
    run_welfare_comparison_efficient_steady_state;
    run_welfare_comparison_inefficient_steady_state;
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Born_Pfeifer_2020
try
    close all; clearvars all; clearvars -global;
    cd('../../Born_Pfeifer_2020');
    run_model_IRF_generation;
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Caldara_et_al_2012
try
    close all; clearvars all; clearvars -global;
    cd('../Caldara_et_al_2012');
    dynare Caldara_et_al_2012
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Chari_et_al_2007
try
    close all; clearvars all; clearvars -global;
    cd('../Chari_et_al_2007');
    dynare Chari_et_al_2007
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Collard_2001
try
    close all; clearvars all; clearvars -global;
    cd('../Collard_2001');
    get_shock_standard_deviation;
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% FV_et_al_2007
try
    close all; clearvars all; clearvars -global;
    cd('../FV_et_al_2007');
    dynare FV_et_al_2007_ABCD
    dynare FV_et_al_2007_ABCD_minreal
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Gali_2008
try
    close all; clearvars all; clearvars -global;
    cd('../Gali_2008');
    dynare Gali_2008_chapter_2
    dynare Gali_2008_chapter_3
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Gali_2010
try
    close all; clearvars all; clearvars -global;
    cd('../Gali_2010');
    dynare Gali_2010
    dynare Gali_2010_calib_target
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Gali_2015
try
    close all; clearvars all; clearvars -global;
    cd('../Gali_2015');
    dynare Gali_2015_chapter_2
    dynare Gali_2015_chapter_3
    dynare Gali_2015_chapter_4
    dynare Gali_2015_chapter_5_discretion
    dynare Gali_2015_chapter_5_commitment
    dynare Gali_2015_chapter_5_discretion_ZLB
    dynare Gali_2015_chapter_5_commitment_ZLB
    dynare Gali_2015_chapter_6
    dynare Gali_2015_chapter_6_4
    dynare Gali_2015_chapter_6_5
    dynare Gali_2015_chapter_7
    dynare Gali_2015_chapter_8
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Gali_Monacelli_2005
try
    close all; clearvars all; clearvars -global;
    cd('../Gali_Monacelli_2005');
    dynare Gali_Monacelli_2005
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% GarciaCicco_et_al_2010
try
    close all; clearvars all; clearvars -global;
    cd('../GarciaCicco_et_al_2010');
    dynare GarciaCicco_et_al_2010
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Ghironi_Melitz_2005
try
    close all; clearvars all; clearvars -global;
    cd('../Ghironi_Melitz_2005');
    dynare Ghironi_Melitz_2005.mod
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Guerrieri_Iacoviello_2015
try
    close all; clearvars all; clearvars -global;
    cd('../Guerrieri_Iacoviello_2015');
    dyn_ver = dynare_version;
    if str2double(dyn_ver(1)) >= 5
        dynare Guerrieri_Iacoviello_2015_rbc
        dynare Guerrieri_Iacoviello_2015_nk
    end
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% HP_filter_missing_data
try
    close all; clearvars all; clearvars -global;
    cd('../HP_filter_missing_data');
    dyn_ver = dynare_version;
    if str2double(dyn_ver(1)) >= 5
        dynare HP_filter_missing_data.mod
    end
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Hansen_1985
try
    close all; clearvars all; clearvars -global;
    cd('../Hansen_1985');
    dynare Hansen_1985.mod
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Ireland_2004
try
    close all; clearvars all; clearvars -global;
    cd('../Ireland_2004');
    dynare Ireland_2004
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Jermann_1998
try
    close all; clearvars all; clearvars -global;
    cd('../Jermann_1998');
    dynare Jermann_1998
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Jermann_Quadrini_2012
try
    close all; clearvars all; clearvars -global;
    cd('../Jermann_Quadrini_2012/Jermann_Quadrini_2012_RBC');
    construct_data
    dynare Jermann_Quadrini_2012_RBC
    cd('../Jermann_Quadrini_2012_NK');
    dynare Jermann_Quadrini_2012_NK
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% McCandless_2008
try
    close all; clearvars all; clearvars -global;
    cd('../../McCandless_2008');
    dynare McCandless_2008_Chapter_9
    dynare McCandless_2008_Chapter_13
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% NK_linear_forward_guidance
try
    close all; clearvars all; clearvars -global;
    cd('../NK_linear_forward_guidance');
    dynare NK_linear_forward_guidance
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% RBC_IRF_matching
try
    close all; clearvars all; clearvars -global;
    cd('../RBC_IRF_matching');
    dynare RBC_IRF_matching
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% RBC_baseline
try
    close all; clearvars all; clearvars -global;
    cd('../RBC_baseline');
    dynare RBC_baseline
    dynare RBC_baseline_first_diff_bayesian
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% RBC_baseline_welfare
try
    close all; clearvars all; clearvars -global;
    cd('../RBC_baseline_welfare');
    dynare RBC_baseline_welfare
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% RBC_capitalstock_shock
try
    close all; clearvars all; clearvars -global;
    cd('../RBC_capitalstock_shock');
    dynare RBC_capitalstock_shock.mod
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% RBC_news_shock_model
try
    close all; clearvars all; clearvars -global;
    cd('../RBC_news_shock_model');
    dynare RBC_news_shock_model
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% RBC_state_dependent_GIRF
try
    close all; clearvars all; clearvars -global;
    cd('../RBC_state_dependent_GIRF');
    dynare RBC_state_dependent_GIRF
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% SGU_2003
try
    close all; clearvars all; clearvars -global;
    cd('../SGU_2003');
    dynare SGU_2003.mod
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% SGU_2004
try
    close all; clearvars all; clearvars -global;
    cd('../SGU_2004');
    dynare SGU_2004
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Sims_2012
try
    close all; clearvars all; clearvars -global;
    cd('../Sims_2012');
    addpath('../FV_et_al_2007'); %ABCD_test.m
    dynare Sims_2012_RBC
    rmpath('../FV_et_al_2007');
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Smets_Wouters_2007
try
    close all; clearvars all; clearvars -global;
    cd('../Smets_Wouters_2007');
    dynare Smets_Wouters_2007
    dynare Smets_Wouters_2007_45
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Solow_model
try
    close all; clearvars all; clearvars -global;
    cd('../Solow_model');
    dynare Solow_SS_transition
    dynare Solow_growth_rate_changes -DTFP_growth=false
    dynare Solow_growth_rate_changes -DTFP_growth=true
    dynare Solow_nonstationary
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Stock_SIR_2020
try
    close all; clearvars all; clearvars -global;
    cd('../Stock_SIR_2020');
    dynare Stock_SIR_2020
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Woodford_2003
try
    close all; clearvars all; clearvars -global;
    cd('../Woodford_2003');
    dynare Woodford_2003_Chapter_7
catch ME
    fid = fopen('error.txt', 'w'); fprintf(fid,'%s',ME.message);fclose(fid);
end

%% Evaluate errors
cd('..')
% Get list of all subfolders.
files = dir(fullfile(pwd, '**', 'error.txt'));
    
% If the file is found, print the folder and error
if ~isempty(files)
    fprintf('\n\n\nFOUND ERRORS:\n***********\n')
    for i = 1 : length(files)
        fprintf('Found error in folder: %s', fullfile(files(i).folder));
        type(fullfile(files(i).folder, files(i).name));
        fprintf('***********\n')
        delete(fullfile(files(i).folder, files(i).name));
    end
    error('There were errors!')
else
    fprintf('NO ERRORS FOUND!\n');
end