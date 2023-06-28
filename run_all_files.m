clearvars all; clearvars -global;

%% Aguiar_Gopinath_2007
close all;
cd('Aguiar_Gopinath_2007');
dynare Aguiar_Gopinath_2007

%% Andreasen_2012
close all; clearvars;
cd('../Andreasen_2012');
dynare Andreasen_2012_rare_disasters

%% Ascari_Sbordone_2014
close all; clearvars;
cd('../Ascari_Sbordone_2014');
dynare Ascari_Sbordone_2014

%% Basu_Bundick_2017
close all; clearvars;
cd('../Basu_Bundick_2017');
dynare Basu_Bundick_2017

%% Born_Pfeifer_2014
close all; clearvars;
cd('../Born_Pfeifer_2014');
dynare Born_Pfeifer_RM_Comment.mod

%% Born_Pfeifer_2018
close all; clearvars;
cd('../Born_Pfeifer_2018/Monetary_Policy_IRFs');
run_IRF_comparison;
cd('../Welfare');
run_welfare_comparison_efficient_steady_state;
run_welfare_comparison_inefficient_steady_state;

%% Born_Pfeifer_2020
close all; clearvars;
cd('../../Born_Pfeifer_2020');
run_model_IRF_generation;

%% Caldara_et_al_2012
close all; clearvars;
cd('../Caldara_et_al_2012');
dynare Caldara_et_al_2012

%% Chari_et_al_2007
close all; clearvars;
cd('../Chari_et_al_2007');
dynare Chari_et_al_2007

%% Collard_2001
close all; clearvars;
cd('../Collard_2001');
get_shock_standard_deviation;

%% FV_et_al_2007
close all; clearvars;
cd('../FV_et_al_2007');
dynare FV_et_al_2007_ABCD
dynare FV_et_al_2007_ABCD_minreal

%% Gali_2008
close all; clearvars;
cd('../Gali_2008');
dynare Gali_2008_chapter_2
dynare Gali_2008_chapter_3

%% Gali_2010
close all; clearvars;
cd('../Gali_2010');
dynare Gali_2010
dynare Gali_2010_calib_target

%% Gali_2015
close all; clearvars;
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

%% Gali_Monacelli_2005
close all; clearvars;
cd('../Gali_Monacelli_2005');
dynare Gali_Monacelli_2005

%% GarciaCicco_et_al_2010
close all; clearvars;
cd('../GarciaCicco_et_al_2010');
dynare GarciaCicco_et_al_2010

%% Ghironi_Melitz_2005
close all; clearvars;
cd('../Ghironi_Melitz_2005');
dynare Ghironi_Melitz_2005.mod

%% Guerrieri_Iacoviello_2015
close all; clearvars;
cd('../Guerrieri_Iacoviello_2015');
dynare Guerrieri_Iacoviello_2015_rbc
dynare Guerrieri_Iacoviello_2015_nk

%% HP_filter_missing_data
close all; clearvars;
cd('../HP_filter_missing_data');
dynare HP_filter_missing_data.mod

%% Hansen_1985
close all; clearvars;
cd('../Hansen_1985');
dynare Hansen_1985.mod

%% Ireland_2004
close all; clearvars;
cd('../Ireland_2004');
dynare Ireland_2004

%% Jermann_1998
close all; clearvars;
cd('../Jermann_1998');
dynare Jermann_1998

%% Jermann_Quadrini_2012
close all; clearvars;
cd('../Jermann_Quadrini_2012/Jermann_Quadrini_2012_RBC');
construct_data
dynare Jermann_Quadrini_2012_RBC
cd('../Jermann_Quadrini_2012_NK');
dynare Jermann_Quadrini_2012_NK

%% McCandless_2008
close all; clearvars;
cd('../../McCandless_2008');
dynare McCandless_2008_Chapter_9
dynare McCandless_2008_Chapter_13

%% NK_linear_forward_guidance
close all; clearvars;
cd('../NK_linear_forward_guidance');
dynare NK_linear_forward_guidance

%% RBC_IRF_matching
close all; clearvars;
cd('../RBC_IRF_matching');
dynare RBC_IRF_matching

%% RBC_baseline
close all; clearvars;
cd('../RBC_baseline');
dynare RBC_baseline
dynare RBC_baseline_first_diff_bayesian

%% RBC_baseline_welfare
close all; clearvars;
cd('../RBC_baseline_welfare');
dynare RBC_baseline_welfare

%% RBC_capitalstock_shock
close all; clearvars;
cd('../RBC_capitalstock_shock');
dynare RBC_capitalstock_shock.mod

%% RBC_news_shock_model
close all; clearvars;
cd('../RBC_news_shock_model');
dynare RBC_news_shock_model

%% RBC_state_dependent_GIRF
close all; clearvars;
cd('../RBC_state_dependent_GIRF');
dynare RBC_state_dependent_GIRF

%% SGU_2003
close all; clearvars;
cd('../SGU_2003');
dynare SGU_2003.mod

%% SGU_2004
close all; clearvars;
cd('../SGU_2004');
dynare SGU_2004

%% Sims_2012
close all; clearvars;
cd('../Sims_2012');
addpath('../FV_et_al_2007'); %ABCD_test.m
dynare Sims_2012_RBC
rmpath('../FV_et_al_2007');

%% Smets_Wouters_2007
close all; clearvars;
cd('../Smets_Wouters_2007');
dynare Smets_Wouters_2007
dynare Smets_Wouters_2007_45

%% Solow_model
close all; clearvars;
cd('../Solow_model');
dynare Solow_SS_transition
dynare Solow_growth_rate_changes
dynare Solow_nonstationary

%% Stock_SIR_2020
close all; clearvars;
cd('../Stock_SIR_2020');
dynare Stock_SIR_2020

%% Woodford_2003
close all; clearvars;
cd('../Woodford_2003');
dynare Woodford_2003_Chapter_7

cd('..')