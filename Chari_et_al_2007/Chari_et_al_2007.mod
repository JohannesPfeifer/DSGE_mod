/*
 * This paper replicates parts of Chari, V.V/Kehoe, Patrick J./McGrattan, Ellen (2007), 
 * Business Cycle Accounting, Econometrica, 75(3), pp. 781–836. 
 * 
 * It demonstrates how to use the linearized benchmark model estimated using Maximum Likelihood to conduct the 
 * Business Cycle Accounting  as is done in the paper for the 1982 recession (the Great Depression exercise relies
 * on non-linear solution techniques). 
 * The mod-file extracts the historical wedges and replicates:
 *      -   Figure 5 - U.S. output and three measured wedges (quarterly, 1979:1–1985:4; normalized to
 *              equal 100 in 1979:1).
 *      -   Table II - Properties of the wedges, 1959:1–2004:3.
 *
 * Notes:
 * - The steady state file is used to translate the estimated parameters into the 
 *      VAR-representation for the exogenous shocks. This makes sure the parameter 
 *      dependencies are correctly taken into account. In particular, the steady state
 *      values for the wedges required for the model equations are translated to the P0
 *      matrix used for the VAR
 * - While the shocks-block accesses parameter values, this does not specify a 
 *      parameter relationship. Rather it sets the respective
 *      entry of the covariance matrix to the value encountered at this stage, i.e. the initial
 *      calibration. Thus, subsequent changes to parameters are not automatically taken into account,
 *      which is the reason a steady state-file is used.
 * - All data are in logs and linearly detrended, while their counterparts in the model are in levels.
 *      The only exception is government spending, which enters in logs in the model.
 * - In the original paper, the treatment of durables services and durables consumption/depreciation in the quarterly 
 *      model is wrong. They enter with a weight of 1% instead of 100% due to omission of the fact 
 *      that the deflator has a base level of 100 (lines 39 and 40 of their mleqtrly/usdata.m). As the durables stock is relatively
 *      acyclical, this is mostly inconsequential for the results. 
 * - The original paper attempts to normalize output to 1 in 1979(1) in the data treatment, but actually 
 *      normalizes it to 1*(1+gamma), which is why log output in 1979Q1, i.e. y(81) is equal 
 *      to log(1.016^(1/4))=0.004
 * - CKM plot the efficiency wedge as the log of TFP, i.e. log(A)=log(exp(z)^(1-theta)) where z is the 
 *      wedge component entering the VAR process
 * - CKM use a penalty function that punishes eigenvalues of the VAR processes for the wedges bigger 
 *      than 0.995 (see footnote a) to Table I, p. 800). Because this cannot be easily done in Dynare 
 *      the penalty is omitted, leading to a somewhat bigger maximum eigenvalue. However, this pushes the 
 *      VAR coefficients closer to the boundary of the stability region. While hardly affecting the point estimates
 *      and conclusion of the paper, it makes the Hessian at the mode not positive definite and renders the standard errors 
 *      invalid.
 *      and the standard errors 
 * - CKM use the linearized model only to extract the investment wedge and the decision rules. All other wedges
 *      are computed based on the original nonlinear model equations. For this purpose, the capital stock is 
 *      initialized at the steady state value in the first period and then iterated forwards. This mod-file also 
 *      shows how to use the Kalman smoother to directly extract the smoothed wedges. As these are based on the 
 *      linearized model, they differ from the ones derived from the nonlinear equations due to Jensen's Inequality.
 *
 * 
 * This implementation was written by Johannes Pfeifer. In
 * case you spot mistakes, email us at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2013-15 Johannes Pfeifer
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <http://www.gnu.org/licenses/>.
 */

%select whether to use the corrected data set where durables enter with
%their full weight (1) or use the original wrong data set (0)
@#define corrected_data=0

%----------------------------------------------------------------
% define variables 
%----------------------------------------------------------------

var y ${y}$ (long_name='output')
    c ${c}$ (long_name='consumption')
    k ${k}$ (long_name='capital')
    x ${x}$ (long_name='investment')
    l ${l}$ (long_name='labor')
    w ${w}$ (long_name='real wage')
    z ${z}$ (long_name='TFP')
    g ${g}$ (long_name='government spending')
    tau_l ${\tau_l}$ (long_name='labor wedge')
    tau_x ${\tau_x}$ (long_name='investment wedge')
    log_labor_wedge ${\log(1-\tau_l)}$ (long_name='log labor wedge')
    log_investment_wedge ${\log\left(\frac{1}{1+\tau_x}\right)}$ (long_name='log investment wedge')
    log_efficiency_wedge ${\log\left(A\right)}$ (long_name='log efficiency wedge')
    ;

varexo eps_z ${\varepsilon^z}$ (long_name='TFP shock') 
    eps_tau_l ${\varepsilon^{\tau_l}}$ (long_name='labor wedge shock')
    eps_tau_x ${\varepsilon^{\tau_x}}$ (long_name='investment wedge shock')
    eps_g ${\varepsilon^g}$ (long_name='government spending shock')
            ;

predetermined_variables k;

%----------------------------------------------------------------
% define parameters
%----------------------------------------------------------------

parameters betta ${\beta}$ (long_name='discount factor')
    beta_hat ${\hat \beta}$ (long_name='discount factor of detrended model')
    theta ${\theta}$ (long_name='capital share')
    psii ${\psi}$ (long_name='labor disutility')
    delta ${\delta}$ (long_name='depreciation rate')
    siggma ${\sigma}$ (long_name='risk aversion')
    gamma_z ${\gamma_z}$ (long_name='technology growth rate')
    gamma_n ${\gamma_n}$ (long_name='population growth rate')
    %Coefficients of the exogenous process; P0 is the mean vector
    P0_z_bar     rho_zz rho_zl rho_zx rho_zg
    P0_tau_l_bar rho_lz rho_ll rho_lx rho_lg
    P0_tau_x_bar rho_xz rho_xl rho_xx rho_xg
    P0_g_bar     rho_gz rho_gl rho_gx rho_gg
    z_bar ${\bar z}$ (long_name='mean TFP')
    tau_l_bar ${\bar {\tau^l}}$ (long_name='mean labor wedge')
    tau_x_bar ${\bar {\tau^x}}$ (long_name='mean investment wedge')
    g_bar$ {\bar g}$ (long_name='mean government spending')    
    sigma_z sigma_tau_l sigma_tau_x sigma_g
    corr_z_tau_l corr_z_tau_x corr_z_g
    corr_tau_l_tau_x corr_tau_l_g
    corr_tau_x_g
    ;

%----------------------------------------------------------------
% set parameter values (as described on p. 798)
%----------------------------------------------------------------
theta   = 0.35; %capital share
psii     = 2.24; %labor disutility
        
% the following are taken from the replication files (mleq.m)
delta   = 1-(1-0.0464)^(1/4); %depreciation rate
betta    = 0.9722^(1/4); %discount factor
beta_hat = 0; %set in steady state; discount factor of detrended model
siggma   = 1; %risk aversion
gamma_n = (1.015)^(1/4)-1; %population growth
gamma_z = (1.016)^(1/4)-1; %TFP growth


%set exogenous processes to values given in paper (Table I)
z_bar = -0.023921304831820;
tau_l_bar = 0.327939686504030;
tau_x_bar = 0.48344057752536;
g_bar = -1.534423337951600;
rho_zz =0.979987375884460; 
rho_zl =-0.013784796530210;
rho_zx = -0.011726790280120;
rho_zg = 0.019238406399280;
rho_lz =-0.032980622374370;
rho_ll = 0.956383107019150;
rho_lx = -0.045084112346010;
rho_lg = 0.056904929264490;
rho_xz = -0.070245267869030;
rho_xl = -0.046005038389430;
rho_xx = 0.896188709044170;
rho_xg = 0.104075465467900;
rho_gz = 0.004810494636130;
rho_gl = -0.008105963085530; 
rho_gx = 0.048839969445590;
rho_gg = 0.971076706106860;       
        
     
verbatim;
%use provided information to compute mean of VAR representation
s_bar=[z_bar;tau_l_bar;tau_x_bar;g_bar];
P=[rho_zz rho_zl rho_zx rho_zg
    rho_lz rho_ll rho_lx rho_lg
    rho_xz rho_xl rho_xx rho_xg
    rho_gz rho_gl rho_gx rho_gg];

P0=(eye(4,4)-P)*s_bar;

%construct correlation matrix from provided Cholesky decomposition of 
Q=[0.011619704018080                   0                   0                   0
   0.001411648230250   0.006440042459250                   0                   0
  -0.010497271810350   0.001031657598730   0.015841640429160                   0
  -0.000575401658550   0.006112475977340   0.014175451027240   0.004583597338150];

V=Q*Q';

Correlation_matrix=diag(sqrt(diag(V)))^(-1)*V*diag(sqrt(diag(V)))^(-1); 
end;

% use computed information on correlation matrix and standard deviations to set parameters used for initialization in the shocks block
sigma_z=sqrt(V(1,1));
sigma_tau_l=sqrt(V(2,2));
sigma_tau_x=sqrt(V(3,3));
sigma_g=sqrt(V(4,4));

corr_z_tau_l=Correlation_matrix(1,2);
corr_z_tau_x=Correlation_matrix(1,3);
corr_z_g=Correlation_matrix(1,4);

corr_tau_l_tau_x=Correlation_matrix(2,3);
corr_tau_l_g=Correlation_matrix(2,4);

corr_tau_x_g=Correlation_matrix(3,4);

P0_z_bar=P0(1);
P0_tau_l_bar=P0(2);
P0_tau_x_bar=P0(3);
P0_g_bar=P0(4);

%----------------------------------------------------------------
% enter model equations
%----------------------------------------------------------------

model; 
    %resource constraint (A.2.1)
    exp(c)+exp(g)+(1+gamma_z)*(1+gamma_n)*exp(k(+1))-(1-delta)*exp(k)=exp(y);
    
    (1+gamma_z)*(1+gamma_n)*exp(k(+1))=(1-delta)*exp(k)+exp(x);
    exp(y)=exp(k)^theta*(exp(z)*exp(l))^(1-theta);
    psii*exp(c)^siggma/(1-exp(l))=(1-tau_l)*exp(w);
    exp(w)=(1-theta)*exp(k)^theta*exp(l)^(-theta)*exp(z)^(1-theta);
    %Euler equation (A.2.3)
    (1+tau_x)*exp(c)^(-siggma)*(1-exp(l))^(psii*(1-siggma))=beta_hat*exp(c(+1))^(-siggma)*(1-exp(l(+1)))^(psii*(1-siggma))*
          (theta*exp(k(+1))^(theta-1)*(exp(z(+1))*exp(l(+1)))^(1-theta)+(1-delta)*(1+tau_x(+1)));
    %observation equations
    % exogenous processes
    z     = P0_z_bar     + rho_zz*z(-1) + rho_zl*tau_l(-1) + rho_zx*tau_x(-1) + rho_zg*g(-1) + eps_z;
    tau_l = P0_tau_l_bar + rho_lz*z(-1) + rho_ll*tau_l(-1) + rho_lx*tau_x(-1) + rho_lg*g(-1) + eps_tau_l;
    tau_x = P0_tau_x_bar + rho_xz*z(-1) + rho_xl*tau_l(-1) + rho_xx*tau_x(-1) + rho_xg*g(-1) + eps_tau_x;
    g     =     P0_g_bar + rho_gz*z(-1) + rho_gl*tau_l(-1) + rho_gx*tau_x(-1) + rho_gg*g(-1) + eps_g;
    log_labor_wedge = log(1-tau_l);
    log_investment_wedge = log(1/(1+tau_x));
    log_efficiency_wedge=log(exp(z)^(1-theta));
end;


%----------------------------------------------------------------
%  set shock variances
%---------------------------------------------------------------

shocks;
var eps_z; stderr sigma_z;
var eps_tau_l; stderr sigma_tau_l;
var eps_tau_x; stderr sigma_tau_x;
var eps_g; stderr sigma_g;

corr eps_z, eps_tau_l =corr_z_tau_l;
corr eps_z, eps_tau_x =corr_z_tau_x;
corr eps_z, eps_g =corr_z_g;

corr eps_tau_l, eps_tau_x =corr_tau_l_tau_x;
corr eps_tau_l, eps_g =corr_tau_l_g;

corr eps_tau_x, eps_g =corr_tau_x_g;
end;


steady;
check;

write_latex_parameter_table;
write_latex_dynamic_model;
collect_latex_files;    

%----------------------------------------------------------------
% generate IRFs
%----------------------------------------------------------------

stoch_simul(order = 1,irf=0,hp_filter=1600);

%%%%%%%%%%%%%%%%%% Replicate original results based on calibrated model and
%%%%%%%%%%%%%%%%%% and linear approximation for tau_x only
        
%%%%%%%% load data

@#if corrected_data==0
    CKM_data=load('Data_CKM_orig');
@#else
    CKM_data=load('Data_CKM_corrected');
@#endif
timeline=1959:0.25:2004.5;
T=length(CKM_data.y);
log_y_t_data = CKM_data.y;
log_x_t_data = CKM_data.x;
log_l_t_data = CKM_data.l;
log_g_t_data = CKM_data.g;

%%%%%%%% get needed steady state values
log_k_ss=oo_.dr.ys(strmatch('k',M_.endo_names,'exact'));
log_x_ss=oo_.dr.ys(strmatch('x',M_.endo_names,'exact'));
log_g_ss=oo_.dr.ys(strmatch('g',M_.endo_names,'exact'));
log_y_ss=oo_.dr.ys(strmatch('y',M_.endo_names,'exact'));
log_c_ss=oo_.dr.ys(strmatch('c',M_.endo_names,'exact'));
log_l_ss=oo_.dr.ys(strmatch('l',M_.endo_names,'exact'));
z_ss=oo_.dr.ys(strmatch('z',M_.endo_names,'exact'));
tau_l_ss=oo_.dr.ys(strmatch('tau_l',M_.endo_names,'exact'));
tau_x_ss=oo_.dr.ys(strmatch('tau_x',M_.endo_names,'exact'));

%%%%%%%%%%% get coefficients of decision rule linking log investment x_t to state variables
        
x_pos_dr_matrices=oo_.dr.inv_order_var(strmatch('x',M_.endo_names,'exact')); %get location of x in rows of dr.gh*
kstate_pos_dr_matrices=find(oo_.dr.state_var==strmatch('k',M_.endo_names,'exact')); %get location of capital k in columns of dr.ghx
x_eps_z_reaction=oo_.dr.ghu(x_pos_dr_matrices,strmatch('eps_z',M_.exo_names,'exact'));
x_eps_tau_l_reaction=oo_.dr.ghu(x_pos_dr_matrices,strmatch('eps_tau_l',M_.exo_names,'exact'));
x_eps_tau_x_reaction=oo_.dr.ghu(x_pos_dr_matrices,strmatch('eps_tau_x',M_.exo_names,'exact'));
x_eps_g_reaction=oo_.dr.ghu(x_pos_dr_matrices,strmatch('eps_g',M_.exo_names,'exact'));
x_k_reaction=oo_.dr.ghx(x_pos_dr_matrices,kstate_pos_dr_matrices);

%%%%%%%% initialize capital stocks
log_k_t=NaN(T+1,1);
K_t=NaN(T+1,1);
log_k_tplus1=NaN(T,1);
K_tplus1=NaN(T,1);
log_k_t(1,1) = log_k_ss;
K_t(1,1)  = exp(log_k_ss);

%%%%%%%% use law of motion for capital to get log capital stock and capital stock
for i=1:T  
    log_k_tplus1(i,1)  = log_k_ss+((1-delta)*(log_k_t(i)-log_k_ss)+exp(log_x_ss)/exp(log_k_ss)*(log_x_t_data(i)-log_x_ss))/((1+gamma_z)*(1+gamma_n));
    log_k_t(i+1,1) = log_k_tplus1(i);
    K_tplus1(i,1)   = ((1-delta)*K_t(i)+exp(log_x_t_data(i)))/((1+gamma_z)*(1+gamma_n));
    K_t(i+1,1)  = K_tplus1(i);
end;
log_k_t      = log_k_t(1:T);
K_t       = K_t(1:T);

%%%%%%%%% compute tau_x_t based on linear solution; for this:
    
% i) compute consumption from resource constraint
log_c_t      = log_c_ss+(exp(log_y_ss)*(log_y_t_data-log_y_ss)-exp(log_x_ss)*(log_x_t_data-log_x_ss)-exp(log_g_ss)*(log_g_t_data-log_g_ss))/exp(log_c_ss);

% ii) compute efficiency wedge from production function
log_z_t      = z_ss+(log_y_t_data-log_y_ss-theta*(log_k_t-log_k_ss))/(1-theta)-(log_l_t_data-log_l_ss);

% iii) compute labor wedge from FOC
tau_l_t    = tau_l_ss+(1-tau_l_ss)*((log_y_t_data-log_y_ss)-(log_c_t-log_c_ss)-1/(1-exp(log_l_ss))*(log_l_t_data-log_l_ss));

% iv) use linear observation equation that relates investment x to states (x=f(tau_l,k,z,g)) to solve for tau_l 
tau_x_t = ((log_x_t_data-log_x_ss)-x_k_reaction*(log_k_t-log_k_ss)-x_eps_z_reaction*(log_z_t-z_ss)-x_eps_tau_l_reaction*(tau_l_t-tau_l_ss)-x_eps_g_reaction*(log_g_t_data-log_g_ss))/x_eps_tau_x_reaction+tau_x_ss;

%%%%%%%%%%%%%% Now that tau_x_t has been recovered, we can use the nonlinear equations to
%%%%%%%%%%%%%% recover the labor wedge and the efficiency wedge:
C_t       = exp(log_y_t_data)-exp(log_x_t_data)-exp(log_g_t_data);
Z_t       = (exp(log_y_t_data)./(K_t.^theta.*exp(log_l_t_data).^(1-theta))).^(1/(1-theta));
Tau_l_t    = 1-psii/(1-theta)* (C_t./exp(log_y_t_data)) .*(exp(log_l_t_data)./(1-exp(log_l_t_data)));

%%%%%%% Define wedges for computations and plotting
        
Output=exp(log_y_t_data-log_y_t_data(81));
Efficiency_wedge  = (Z_t/Z_t(81)).^(1-theta);
Labor_wedge  = (1-Tau_l_t)/(1-Tau_l_t(81));
Investment_wedge  = (1+tau_x_t(81))./(1+tau_x_t);
Government_wedge  = exp(log_g_t_data-log_g_t_data(81));

figure('Name','Figure 5 (Replication)')
orient portrait
plot(timeline,100*Output,'b-',...
        timeline,100*Efficiency_wedge,'--', ...
        timeline,100*Labor_wedge,':', ...
        timeline,100*Investment_wedge,'-.','LineWidth',1.5)
hold on
plot(timeline,100*ones(T,1),'k:')
title('Figure 5. U.S. output and three measured Wedges (quarterly, 1979:1-1985:4), 1979:1=100')
legend('Data','Efficiency Wedge','Labor Wedge','Investment Wedge','Location','NorthWest')
axis([1979,1986,90,110])

%%%%%% get moments of wedges
        
[~, y_cyc]=sample_hp_filter(log(Output),1600);
[~, log_efficiency_wedge_cyc]=sample_hp_filter(log(Efficiency_wedge),1600);
[~, log_labor_wedge_cyc]=sample_hp_filter(log(Labor_wedge),1600);
[~, log_investment_wedge_cyc]=sample_hp_filter(log(Investment_wedge),1600);
[~, log_government_wedge_cyc]=sample_hp_filter(log(Government_wedge),1600);

relative_volatility(1,1)=std(log_efficiency_wedge_cyc)/std(y_cyc);
relative_volatility(2,1)=std(log_labor_wedge_cyc)/std(y_cyc);
relative_volatility(3,1)=std(log_investment_wedge_cyc)/std(y_cyc);
relative_volatility(4,1)=std(log_government_wedge_cyc)/std(y_cyc);

wedge_data_mat=[log_efficiency_wedge_cyc,log_labor_wedge_cyc,log_investment_wedge_cyc,log_government_wedge_cyc];
cross_correlations=NaN(4,5);
cross_correlations(:,1)=corr(y_cyc(1+2:end),wedge_data_mat(1:end-2,:));
cross_correlations(:,2)=corr(y_cyc(1+1:end),wedge_data_mat(1:end-1,:));
cross_correlations(:,3)=corr(y_cyc,wedge_data_mat);
cross_correlations(:,4)=corr(y_cyc(1:end-1),wedge_data_mat(1+1:end,:));
cross_correlations(:,5)=corr(y_cyc(1:end-2),wedge_data_mat(1+2:end,:));

wedge_names={'Efficiency_wedge','Labor_wedge','Investment_wedge','Government_wedge'};
fprintf('\nTable II:PROPERTIES OF THE WEDGES, 1959:1–2004:3\n')
fprintf('Panel a: Summary statistics\n')
fprintf('%-30s \t %7s \t %7s \t %7s \t %7s \t %7s \t %7s \n','Wedges','Rel. std','-2 ','-1 ','0 ','1 ','2 ')

for ii=1:4
    fprintf('%-30s \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \n',wedge_names{ii},relative_volatility(ii),cross_correlations(ii,:));
end

cross_correlations=NaN(6,4);
wedge_names={'Efficiency Wedge','Labor Wedge','Investment Wedge','Government Wedge'};
iter=1;
for ii=1:3
    for jj=ii+1:4
        row_header{iter,1}=[wedge_names{ii},', ',wedge_names{jj}];
        cross_correlations(iter,5)=corr(wedge_data_mat(1+2:end,ii),wedge_data_mat(1:end-2,jj));
        cross_correlations(iter,4)=corr(wedge_data_mat(1+1:end,ii),wedge_data_mat(1:end-1,jj));
        cross_correlations(iter,3)=corr(wedge_data_mat(:,ii),wedge_data_mat(:,jj));
        cross_correlations(iter,2)=corr(wedge_data_mat(1:end-1,ii),wedge_data_mat(1+1:end,jj));
        cross_correlations(iter,1)=corr(wedge_data_mat(1:end-2,ii),wedge_data_mat(1+2:end,jj));
        iter=iter+1;
    end
end

fprintf('Panel b: Cross-correlations\n')
for ii=1:length(row_header)
    fprintf('%-40s \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \n',row_header{ii},cross_correlations(ii,:));
end

%%%%%%%%%%%%%%%%%%%%%%%% Create same figure with smoothed wedges from Kalman 
%%%%%%%%%%%%%%%%%%%%%%%% smoother, i.e. based on purely linear model
        
log_efficiency_wedge_pos=strmatch('log_efficiency_wedge',M_.endo_names,'exact');
log_g_pos=strmatch('g',M_.endo_names,'exact');
log_labor_wedge_pos=strmatch('log_labor_wedge',M_.endo_names,'exact');
log_investment_wedge_pos=strmatch('log_investment_wedge',M_.endo_names,'exact');
log_y_pos=strmatch('y',M_.endo_names,'exact');

// Theoretical standard deviations
// relative_standard_deviations=NaN(4,1);
// relative_standard_deviations(1,1)=sqrt(oo_.var(log_efficiency_wedge_pos,log_efficiency_wedge_pos))/sqrt(oo_.var(log_y_pos,log_y_pos));
// relative_standard_deviations(2,1)=sqrt(oo_.var(log_labor_wedge_pos,log_labor_wedge_pos))/sqrt(oo_.var(log_y_pos,log_y_pos));
// relative_standard_deviations(3,1)=sqrt(oo_.var(log_investment_wedge_pos,log_investment_wedge_pos))/sqrt(oo_.var(log_y_pos,log_y_pos));
// relative_standard_deviations(4,1)=sqrt(oo_.var(log_g_pos,log_g_pos))/sqrt(oo_.var(log_y_pos,log_y_pos));

varobs y x l g;

%%%%%%%%%%%%%%%%%% run smoother to extract unobserved states, based on linear model
@#if corrected_data==0
    calib_smoother(datafile=Data_CKM_orig);
@#else
    calib_smoother(datafile=Data_CKM_corrected);
@#endif

normalized_log_output=oo_.SmoothedVariables.y-oo_.SmoothedVariables.y(81)+1;
[~, y_cyc]=sample_hp_filter(normalized_log_output,1600);

normalized_log_efficiency_wedge=oo_.SmoothedVariables.log_efficiency_wedge-oo_.SmoothedVariables.log_efficiency_wedge(81)+1;
[~, log_efficiency_wedge_cyc]=sample_hp_filter(normalized_log_efficiency_wedge,1600);

normalized_log_labor_wedge=oo_.SmoothedVariables.log_labor_wedge-oo_.SmoothedVariables.log_labor_wedge(81)+1;
[~, log_labor_wedge_cyc]=sample_hp_filter(normalized_log_labor_wedge,1600);

normalized_log_investment_wedge=oo_.SmoothedVariables.log_investment_wedge-oo_.SmoothedVariables.log_investment_wedge(81)+1;
[~, log_investment_wedge_cyc]=sample_hp_filter(normalized_log_investment_wedge,1600);

normalized_log_government_wedge=oo_.SmoothedVariables.g-oo_.SmoothedVariables.g(81)+1;
[~, log_government_wedge_cyc]=sample_hp_filter(normalized_log_government_wedge,1600);

relative_volatility(1,1)=std(log_efficiency_wedge_cyc)/std(y_cyc);
relative_volatility(2,1)=std(log_labor_wedge_cyc)/std(y_cyc);
relative_volatility(3,1)=std(log_investment_wedge_cyc)/std(y_cyc);
relative_volatility(4,1)=std(log_government_wedge_cyc)/std(y_cyc);


timeline=1959:0.25:2004.5;
condition=(timeline>=1979 & timeline<1986);

figure('Name','Figure 5 (linearized model)')
efficiency_wedge=exp(oo_.SmoothedVariables.log_efficiency_wedge(condition));
plot(timeline(condition),efficiency_wedge./efficiency_wedge(1)*100,'--','LineWidth',1.5)
hold on
labor_wedge=1-oo_.SmoothedVariables.tau_l(condition);
plot(timeline(condition),labor_wedge./labor_wedge(1)*100,':','LineWidth',1.5)
investment_wedge=1./(1+oo_.SmoothedVariables.tau_x(condition));
plot(timeline(condition),investment_wedge./investment_wedge(1)*100,'-.','LineWidth',1.5)

CKM_data=load('Data_CKM_orig');
output=CKM_data.y(condition);
plot(timeline(condition),exp(output-output(1))*100,'b-','LineWidth',1.5)
plot(timeline(condition),100*ones(length(investment_wedge),1),'-')
legend('Efficiency Wedge','Labor Wege','Investment Wedge','Output')
axis tight

estimated_params;
P0_z_bar,  P0_z_bar;
P0_tau_l_bar,  P0_tau_l_bar;
P0_tau_x_bar,  P0_tau_x_bar;
P0_g_bar,  P0_g_bar;
rho_zz, 0.980; 
rho_zl, -0.0138;
rho_zx,  -0.0117;
rho_zg,  0.0192;
rho_lz, -0.0330;
rho_ll,  0.956;
rho_lx,  -0.0451;
rho_lg,  0.0569;
rho_xz,  -0.0702;
rho_xl,  -0.0460;
rho_xx,  0.896;
rho_xg,  0.104;
rho_gz,  0.00481;
rho_gl,  -0.00811; 
rho_gx,  0.0488;
rho_gg,  0.971;
stderr  eps_z, sigma_z;
stderr  eps_tau_l, sigma_tau_l;
stderr  eps_tau_x, sigma_tau_x;
stderr  eps_g, sigma_g;

corr  eps_z, eps_tau_l, corr_z_tau_l;
corr  eps_z, eps_tau_x, corr_z_tau_x;
corr  eps_z, eps_g, corr_z_g;

corr  eps_tau_l, eps_tau_x, corr_tau_l_tau_x;
corr  eps_tau_l, eps_g, corr_tau_l_g;

corr  eps_tau_x, eps_g, corr_tau_x_g;
end;

@#if corrected_data==0
    estimation(datafile=Data_CKM_orig,mode_compute=4,mode_check,smoother,consider_all_endogenous,silent_optimizer);
@#else
    estimation(datafile=Data_CKM_corrected,mode_compute=4,mode_check,smoother,consider_all_endogenous,silent_optimizer);
@#endif

verbatim;
s_bar=[M_.params(strmatch('z_bar',deblank(M_.param_names),'exact'));M_.params(strmatch('tau_l_bar',deblank(M_.param_names),'exact'));M_.params(strmatch('tau_x_bar',deblank(M_.param_names),'exact'));M_.params(strmatch('g_bar',deblank(M_.param_names),'exact'))];
P=[M_.params(strmatch('rho_zz',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_zl',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_zx',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_zg',deblank(M_.param_names),'exact'))
    M_.params(strmatch('rho_lz',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_ll',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_lz',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_lg',deblank(M_.param_names),'exact'))
    M_.params(strmatch('rho_xz',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_xl',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_xx',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_xg',deblank(M_.param_names),'exact'))
    M_.params(strmatch('rho_gz',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_gl',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_gx',deblank(M_.param_names),'exact')) M_.params(strmatch('rho_gg',deblank(M_.param_names),'exact'))];
P0=(eye(4,4)-P)*s_bar;
Q=chol(M_.Sigma_e,'lower');
end;

%%%%%%%% get needed steady state values
log_k_ss=oo_.dr.ys(strmatch('k',M_.endo_names,'exact'));
log_x_ss=oo_.dr.ys(strmatch('x',M_.endo_names,'exact'));
log_g_ss=oo_.dr.ys(strmatch('g',M_.endo_names,'exact'));
log_y_ss=oo_.dr.ys(strmatch('y',M_.endo_names,'exact'));
log_c_ss=oo_.dr.ys(strmatch('c',M_.endo_names,'exact'));
log_l_ss=oo_.dr.ys(strmatch('l',M_.endo_names,'exact'));
z_ss=oo_.dr.ys(strmatch('z',M_.endo_names,'exact'));
tau_l_ss=oo_.dr.ys(strmatch('tau_l',M_.endo_names,'exact'));
tau_x_ss=oo_.dr.ys(strmatch('tau_x',M_.endo_names,'exact'));

%%%%%%%%%%% get coefficients of decision rule linking log investment x_t to state variables
        
x_pos_dr_matrices=oo_.dr.inv_order_var(strmatch('x',M_.endo_names,'exact')); %get location of x in rows of dr.gh*
kstate_pos_dr_matrices=find(oo_.dr.state_var==strmatch('k',M_.endo_names,'exact')); %get location of capital k in columns of dr.ghx
x_eps_z_reaction=oo_.dr.ghu(x_pos_dr_matrices,strmatch('eps_z',M_.exo_names,'exact'));
x_eps_tau_l_reaction=oo_.dr.ghu(x_pos_dr_matrices,strmatch('eps_tau_l',M_.exo_names,'exact'));
x_eps_tau_x_reaction=oo_.dr.ghu(x_pos_dr_matrices,strmatch('eps_tau_x',M_.exo_names,'exact'));
x_eps_g_reaction=oo_.dr.ghu(x_pos_dr_matrices,strmatch('eps_g',M_.exo_names,'exact'));
x_k_reaction=oo_.dr.ghx(x_pos_dr_matrices,kstate_pos_dr_matrices);

%%%%%%%% initialize capital stocks
log_k_t=NaN(T+1,1);
K_t=NaN(T+1,1);
log_k_tplus1=NaN(T,1);
K_tplus1=NaN(T,1);
log_k_t(1,1) = log_k_ss;
K_t(1,1)  = exp(log_k_ss);

%%%%%%%% use law of motion for capital to get log capital stock and capital stock
for i=1:T  
    log_k_tplus1(i,1)  = log_k_ss+((1-delta)*(log_k_t(i)-log_k_ss)+exp(log_x_ss)/exp(log_k_ss)*(log_x_t_data(i)-log_x_ss))/((1+gamma_z)*(1+gamma_n));
    log_k_t(i+1,1) = log_k_tplus1(i);
    K_tplus1(i,1)   = ((1-delta)*K_t(i)+exp(log_x_t_data(i)))/((1+gamma_z)*(1+gamma_n));
    K_t(i+1,1)  = K_tplus1(i);
end;
log_k_t      = log_k_t(1:T);
K_t       = K_t(1:T);

%%%%%%%%% compute tau_x_t based on linear solution; for this:
    
% i) compute consumption from resource constraint
log_c_t      = log_c_ss+(exp(log_y_ss)*(log_y_t_data-log_y_ss)-exp(log_x_ss)*(log_x_t_data-log_x_ss)-exp(log_g_ss)*(log_g_t_data-log_g_ss))/exp(log_c_ss);

% ii) compute efficiency wedge from production function
log_z_t      = z_ss+(log_y_t_data-log_y_ss-theta*(log_k_t-log_k_ss))/(1-theta)-(log_l_t_data-log_l_ss);

% iii) compute labor wedge from FOC
tau_l_t    = tau_l_ss+(1-tau_l_ss)*((log_y_t_data-log_y_ss)-(log_c_t-log_c_ss)-1/(1-exp(log_l_ss))*(log_l_t_data-log_l_ss));

% iv) use linear observation equation that relates investment x to states (x=f(tau_l,k,z,g)) to solve for tau_l 
tau_x_t = ((log_x_t_data-log_x_ss)-x_k_reaction*(log_k_t-log_k_ss)-x_eps_z_reaction*(log_z_t-z_ss)-x_eps_tau_l_reaction*(tau_l_t-tau_l_ss)-x_eps_g_reaction*(log_g_t_data-log_g_ss))/x_eps_tau_x_reaction+tau_x_ss;

%%%%%%%%%%%%%% Now that tau_x_t has been recovered, we can use the nonlinear equations to
%%%%%%%%%%%%%% recover the labor wedge and the efficiency wedge:
C_t       = exp(log_y_t_data)-exp(log_x_t_data)-exp(log_g_t_data);
Z_t       = (exp(log_y_t_data)./(K_t.^theta.*exp(log_l_t_data).^(1-theta))).^(1/(1-theta));
Tau_l_t    = 1-psii/(1-theta)* (C_t./exp(log_y_t_data)) .*(exp(log_l_t_data)./(1-exp(log_l_t_data)));

%%%%%%% Define wedges for computations and plotting
        
Output=exp(log_y_t_data-log_y_t_data(81));
Efficiency_wedge  = (Z_t/Z_t(81)).^(1-theta);
Labor_wedge  = (1-Tau_l_t)/(1-Tau_l_t(81));
Investment_wedge  = (1+tau_x_t(81))./(1+tau_x_t);
Government_wedge  = exp(log_g_t_data-log_g_t_data(81));

figure('Name','Figure 5 (Replication)')
orient portrait
plot(timeline,100*Output,'b-',...
        timeline,100*Efficiency_wedge,'--', ...
        timeline,100*Labor_wedge,':', ...
        timeline,100*Investment_wedge,'-.','LineWidth',1.5)
hold on
plot(timeline,100*ones(T,1),'k:')
title('Figure 5. U.S. output and three measured Wedges (quarterly, 1979:1-1985:4), 1979:1=100')
legend('Data','Efficiency Wedge','Labor Wedge','Investment Wedge','Location','NorthWest')
axis([1979,1986,90,110])

%%%%%% get moments of wedges
        
[~, y_cyc]=sample_hp_filter(log(Output),1600);
[~, log_efficiency_wedge_cyc]=sample_hp_filter(log(Efficiency_wedge),1600);
[~, log_labor_wedge_cyc]=sample_hp_filter(log(Labor_wedge),1600);
[~, log_investment_wedge_cyc]=sample_hp_filter(log(Investment_wedge),1600);
[~, log_government_wedge_cyc]=sample_hp_filter(log(Government_wedge),1600);

relative_volatility(1,1)=std(log_efficiency_wedge_cyc)/std(y_cyc);
relative_volatility(2,1)=std(log_labor_wedge_cyc)/std(y_cyc);
relative_volatility(3,1)=std(log_investment_wedge_cyc)/std(y_cyc);
relative_volatility(4,1)=std(log_government_wedge_cyc)/std(y_cyc);

wedge_data_mat=[log_efficiency_wedge_cyc,log_labor_wedge_cyc,log_investment_wedge_cyc,log_government_wedge_cyc];
cross_correlations=NaN(4,5);
cross_correlations(:,1)=corr(y_cyc(1+2:end),wedge_data_mat(1:end-2,:));
cross_correlations(:,2)=corr(y_cyc(1+1:end),wedge_data_mat(1:end-1,:));
cross_correlations(:,3)=corr(y_cyc,wedge_data_mat);
cross_correlations(:,4)=corr(y_cyc(1:end-1),wedge_data_mat(1+1:end,:));
cross_correlations(:,5)=corr(y_cyc(1:end-2),wedge_data_mat(1+2:end,:));

wedge_names={'Efficiency_wedge','Labor_wedge','Investment_wedge','Government_wedge'};
fprintf('\nTable II:PROPERTIES OF THE WEDGES, 1959:1–2004:3\n')
fprintf('Panel a: Summary statistics\n')
fprintf('%-30s \t %7s \t %7s \t %7s \t %7s \t %7s \t %7s \n','Wedges','Rel. std','-2 ','-1 ','0 ','1 ','2 ')

for ii=1:4
    fprintf('%-30s \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \n',wedge_names{ii},relative_volatility(ii),cross_correlations(ii,:));
end

cross_correlations=NaN(6,4);
wedge_names={'Efficiency Wedge','Labor Wedge','Investment Wedge','Government Wedge'};
iter=1;
for ii=1:3
    for jj=ii+1:4
        row_header{iter,1}=[wedge_names{ii},', ',wedge_names{jj}];
        cross_correlations(iter,5)=corr(wedge_data_mat(1+2:end,ii),wedge_data_mat(1:end-2,jj));
        cross_correlations(iter,4)=corr(wedge_data_mat(1+1:end,ii),wedge_data_mat(1:end-1,jj));
        cross_correlations(iter,3)=corr(wedge_data_mat(:,ii),wedge_data_mat(:,jj));
        cross_correlations(iter,2)=corr(wedge_data_mat(1:end-1,ii),wedge_data_mat(1+1:end,jj));
        cross_correlations(iter,1)=corr(wedge_data_mat(1:end-2,ii),wedge_data_mat(1+2:end,jj));
        iter=iter+1;
    end
end

fprintf('Panel b: Cross-correlations\n')
for ii=1:length(row_header)
    fprintf('%-40s \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \t % 4.3f \n',row_header{ii},cross_correlations(ii,:));
end

@#if corrected_data==0
    save CKM_Results_original_data
@#else
    save CKM_Results_corrected_data
@#endif