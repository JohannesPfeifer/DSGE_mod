%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'Jermann_1998';
M_.dynare_version = '4.6-unstable-a51fef393a44616e399041870983ad450945cab3';
oo_.dynare_version = '4.6-unstable-a51fef393a44616e399041870983ad450945cab3';
options_.dynare_version = '4.6-unstable-a51fef393a44616e399041870983ad450945cab3';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('Jermann_1998.log');
M_.exo_names = 'e';
M_.exo_names_tex = '{\beta^*}';
M_.exo_names_long = 'productivity shock';
M_.endo_names = 'c';
M_.endo_names_tex = 'C';
M_.endo_names_long = 'consumption';
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'K');
M_.endo_names_long = char(M_.endo_names_long, 'capital stock');
M_.endo_names = char(M_.endo_names, 'invest');
M_.endo_names_tex = char(M_.endo_names_tex, 'I');
M_.endo_names_long = char(M_.endo_names_long, 'investment');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'productivity');
M_.endo_names = char(M_.endo_names, 'lambda');
M_.endo_names_tex = char(M_.endo_names_tex, '\lambda');
M_.endo_names_long = char(M_.endo_names_long, 'marginal utility');
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'Tobins marginal q');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'W');
M_.endo_names_long = char(M_.endo_names_long, 'real wage');
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'dividends');
M_.endo_names = char(M_.endo_names, 'V_k');
M_.endo_names_tex = char(M_.endo_names_tex, '{V^k}');
M_.endo_names_long = char(M_.endo_names_long, 'price of capital stock');
M_.endo_names = char(M_.endo_names, 'V_b');
M_.endo_names_tex = char(M_.endo_names_tex, '{V^b}');
M_.endo_names_long = char(M_.endo_names_long, 'price of a consol');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y');
M_.endo_names_long = char(M_.endo_names_long, 'output');
M_.endo_names = char(M_.endo_names, 'r_k');
M_.endo_names_tex = char(M_.endo_names_tex, '{r^K}');
M_.endo_names_long = char(M_.endo_names_long, 'return to capital');
M_.endo_names = char(M_.endo_names, 'r_f');
M_.endo_names_tex = char(M_.endo_names_tex, '{r^f}');
M_.endo_names_long = char(M_.endo_names_long, 'risk-free interest rate');
M_.endo_names = char(M_.endo_names, 'c_growth');
M_.endo_names_tex = char(M_.endo_names_tex, '{\Delta C}');
M_.endo_names_long = char(M_.endo_names_long, 'growth rate of consumption');
M_.endo_names = char(M_.endo_names, 'y_growth');
M_.endo_names_tex = char(M_.endo_names_tex, '{\Delta Y}');
M_.endo_names_long = char(M_.endo_names_long, 'growth rate of output');
M_.endo_names = char(M_.endo_names, 'i_growth');
M_.endo_names_tex = char(M_.endo_names_tex, '{\Delta I}');
M_.endo_names_long = char(M_.endo_names_long, 'growth rate of investment');
M_.endo_names = char(M_.endo_names, 'log_d');
M_.endo_names_tex = char(M_.endo_names_tex, '{\log(d)}');
M_.endo_names_long = char(M_.endo_names_long, 'log dividend');
M_.endo_names = char(M_.endo_names, 'log_V_k');
M_.endo_names_tex = char(M_.endo_names_tex, '{\log(V^k)}');
M_.endo_names_long = char(M_.endo_names_long, 'log price of capital stock');
M_.endo_names = char(M_.endo_names, 'log_V_b');
M_.endo_names_tex = char(M_.endo_names_tex, '{\log(V^b)}');
M_.endo_names_long = char(M_.endo_names_long, 'log price of a consol');
M_.endo_names = char(M_.endo_names, 'log_r_f');
M_.endo_names_tex = char(M_.endo_names_tex, '{\log(r^f)}');
M_.endo_names_long = char(M_.endo_names_long, 'log risk-free rate');
M_.endo_names = char(M_.endo_names, 'log_r_k');
M_.endo_names_tex = char(M_.endo_names_tex, '{\log(r^k)}');
M_.endo_names_long = char(M_.endo_names_long, 'log return to capital');
M_.endo_names = char(M_.endo_names, 'log_lambda');
M_.endo_names_tex = char(M_.endo_names_tex, '{\log(\lambda)}');
M_.endo_names_long = char(M_.endo_names_long, 'log marginal utility');
M_.endo_names = char(M_.endo_names, 'SDF');
M_.endo_names_tex = char(M_.endo_names_tex, '{\beta\frac{U_{C,t+1}}{U_{C,t}}}');
M_.endo_names_long = char(M_.endo_names_long, 'stochastic discount factor');
M_.endo_names = char(M_.endo_names, 'rf_ann');
M_.endo_names_tex = char(M_.endo_names_tex, '{rf^{ann,net}}');
M_.endo_names_long = char(M_.endo_names_long, 'annualized net risk free rate');
M_.endo_names = char(M_.endo_names, 'rk_ann');
M_.endo_names_tex = char(M_.endo_names_tex, '{rk^{ann,net}}');
M_.endo_names_long = char(M_.endo_names_long, 'annualized net return to capital');
M_.endo_names = char(M_.endo_names, 'rp_ann');
M_.endo_names_tex = char(M_.endo_names_tex, '{ERP^{ann}}');
M_.endo_names_long = char(M_.endo_names_long, 'annualized equity premium');
M_.endo_names = char(M_.endo_names, 'equity_premium');
M_.endo_names_tex = char(M_.endo_names_tex, '{ERP}');
M_.endo_names_long = char(M_.endo_names_long, 'equity premium');
M_.endo_partitions = struct();
M_.param_names = 'betastar';
M_.param_names_tex = '{\beta^*}';
M_.param_names_long = 'discount factor';
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, '{\delta}');
M_.param_names_long = char(M_.param_names_long, 'depreciation rate');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, '{\alpha}');
M_.param_names_long = char(M_.param_names_long, 'capital share');
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, '{\sigma}');
M_.param_names_long = char(M_.param_names_long, 'standard deviation productivity shock');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, '{\rho}');
M_.param_names_long = char(M_.param_names_long, 'autocorrelation productivity');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, '{\gamma}');
M_.param_names_long = char(M_.param_names_long, 'long-run growth rate');
M_.param_names = char(M_.param_names, 'tau');
M_.param_names_tex = char(M_.param_names_tex, '{\tau}');
M_.param_names_long = char(M_.param_names_long, 'risk aversion');
M_.param_names = char(M_.param_names, 'h');
M_.param_names_tex = char(M_.param_names_tex, '{h}');
M_.param_names_long = char(M_.param_names_long, 'habit parameter (originally called alpha in the paper)');
M_.param_names = char(M_.param_names, 'xi');
M_.param_names_tex = char(M_.param_names_tex, '{\xi}');
M_.param_names_long = char(M_.param_names_long, 'elasticity of investment-capital ratio wrt to q');
M_.param_names = char(M_.param_names, 'const');
M_.param_names_tex = char(M_.param_names_tex, '{c}');
M_.param_names_long = char(M_.param_names_long, 'constant in in investment adjustment costs');
M_.param_names = char(M_.param_names, 'a');
M_.param_names_tex = char(M_.param_names_tex, '{a}');
M_.param_names_long = char(M_.param_names_long, 'Parameter a in investment adjustment costs');
M_.param_names = char(M_.param_names, 'b');
M_.param_names_tex = char(M_.param_names_tex, '{b}');
M_.param_names_long = char(M_.param_names_long, 'Parameter b in investment adjustment costs');
M_.param_names = char(M_.param_names, 'i_k');
M_.param_names_tex = char(M_.param_names_tex, '{\frac{I}{K}}');
M_.param_names_long = char(M_.param_names_long, 'steady state investment to capital ratio');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 27;
M_.param_nbr = 13;
M_.orig_endo_nbr = 27;
M_.aux_vars = [];
M_.predetermined_variables = [ 2 ];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 0;
erase_compiled_function('Jermann_1998_static');
erase_compiled_function('Jermann_1998_dynamic');
M_.orig_eq_nbr = 27;
M_.eq_nbr = 27;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.max_endo_lag_orig = 1;
M_.max_endo_lead_orig = 1;
M_.max_exo_lag_orig = 0;
M_.max_exo_lead_orig = 0;
M_.max_exo_det_lag_orig = 0;
M_.max_exo_det_lead_orig = 0;
M_.max_lag_orig = 1;
M_.max_lead_orig = 1;
M_.lead_lag_incidence = [
 1 6 33;
 2 7 0;
 3 8 34;
 4 9 35;
 0 10 36;
 0 11 37;
 0 12 0;
 0 13 38;
 0 14 39;
 0 15 40;
 5 16 0;
 0 17 0;
 0 18 0;
 0 19 0;
 0 20 0;
 0 21 0;
 0 22 0;
 0 23 0;
 0 24 0;
 0 25 0;
 0 26 0;
 0 27 0;
 0 28 0;
 0 29 0;
 0 30 0;
 0 31 0;
 0 32 0;]';
M_.nstatic = 17;
M_.nfwrd   = 5;
M_.npred   = 2;
M_.nboth   = 3;
M_.nsfwrd   = 8;
M_.nspred   = 5;
M_.ndynamic   = 10;
M_.equations_tags = {
  1 , 'name' , '1. Marginal utility' ;
  2 , 'name' , '2. Euler equation stocks' ;
  3 , 'name' , 'Definition dividends' ;
  4 , 'name' , 'Real wage' ;
  5 , 'name' , 'Resource constraint' ;
  6 , 'name' , 'LOM capital' ;
  7 , 'name' , 'FOC capital' ;
  8 , 'name' , 'FOC investment' ;
  9 , 'name' , 'LOM technology' ;
  10 , 'name' , 'Production function' ;
  11 , 'name' , 'Return to capital' ;
  12 , 'name' , 'Definition stochastic discount factor' ;
  13 , 'name' , 'Definition risk-free interest rate' ;
  14 , 'name' , 'log risk-free interest rate' ;
  15 , 'name' , 'log return to capital' ;
  16 , 'name' , 'log return to capital' ;
  17 , 'name' , 'Annualized net risk-free rate' ;
  18 , 'name' , 'Annualized net return to capital' ;
  19 , 'name' , 'log return to capital' ;
  20 , 'name' , 'Pricing equation for a consol' ;
  21 , 'name' , 'Definition growth rate of output' ;
  22 , 'name' , 'Definition growth rate of investment' ;
  23 , 'name' , 'Definition growth rate of consumption' ;
  24 , 'name' , 'Definition log price of capital' ;
  25 , 'name' , 'Definition log price of consol' ;
  26 , 'name' , 'Definition log dividend' ;
  27 , 'name' , 'Definition log marginal utility' ;
};
M_.static_and_dynamic_models_differ = 0;
M_.state_var = [1 2 3 4 11 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(27, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(13, 1);
M_.NNZDerivatives = [85; 102; -1];
M_.params( 6 ) = 1.005;
gamma = M_.params( 6 );
M_.params( 3 ) = 0.36;
alpha = M_.params( 3 );
M_.params( 2 ) = 0.025;
delta = M_.params( 2 );
M_.params( 7 ) = 5;
tau = M_.params( 7 );
M_.params( 8 ) = 0.82;
h = M_.params( 8 );
M_.params( 1 ) = M_.params(6)/1.011138;
betastar = M_.params( 1 );
M_.params( 4 ) = 0.0064/(1-M_.params(3));
sigma = M_.params( 4 );
M_.params( 5 ) = 0.99;
rho = M_.params( 5 );
M_.params( 9 ) = .23;
xi = M_.params( 9 );
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(4)^2;
steady;
options_.irf = 0;
options_.order = 1;
var_list_ = char('y_growth','c_growth','i_growth');
info = stoch_simul(var_list_);
options_.irf = 0;
options_.order = 2;
options_.periods = 50000;
var_list_ = char('V_k','d','r_f','r_k');
info = stoch_simul(var_list_);
E_r_f=mean(exp(log_r_f)-1)*400
R=gamma*(exp(log_V_k(2:end))+exp(log_d(2:end)))./exp(log_V_k(1:end-1));
E_r_k=(mean(R)-1)*400
R_b=gamma*(exp(log_V_b(2:end))+1/betastar-1)./exp(log_V_b(1:end-1));
E_r_b=(mean(R_b)-1)*400
E_r_k-E_r_b
save('Jermann_1998_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('Jermann_1998_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('Jermann_1998_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('Jermann_1998_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('Jermann_1998_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('Jermann_1998_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('Jermann_1998_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
disp('Note: 2 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
