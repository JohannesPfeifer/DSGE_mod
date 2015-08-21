/*
 * This file replicates the model studied in:
 * Ascari, Guido and Sbordone, Argia M. (2014): "The Macroeconomics of Trend Inflation",
 * Journal of Economic Literature, 52(3), pp. 679-739.
 * 
 * It provides a replication of the main results of the original paper  
 * in Section 3 (the New Keynesian model with trend inflation). It replicates the Figures:
 * - Figure 7: The Cost of Price dispersion
 * - Figure 8: Trend Inflation and Steady State Variables
 * - Figure 11: The Determinacy Region and Trend Inflation
 * - Figure 13: Impulse Response Functions to a 1 Percent Positive Technology Shock 
 * - Figure 14: Impulse Response Functions to a 1 Percent Positive Monetary Policy Shock 
 *
 * Moreover, it replicates the business cycle moments reported on p. 717. 
 *
 * This mod-file shows how to access steady state variables in order to plot steady state
 * dependences on parameters.
 * It also shows how to manually do a stability mapping be iterating over a grid on the parameter space.
 * 
 * Notes:
 * - The mod-file requires Dynare > 4.4.3, i.e. currently the unstable version. Otherwise, you will get error like 
 *      "ERROR: in the 'steady_state' block, variable 'pi' is undefined in the declaration of variable 'pi'"
 * - The results from the nonlinear model have been cross-checked with the linearized 
 *      version presented in the paper and the replication files provided by the authors
 * - The annual inflation target is disaggregated to quarterly figures using the geometric mean. Simply dividing by
 *      4 results in small numerical differences
 * - Following the approach of the published replication files by the authors, the labor disutility parameter, which
 *      is unspecified in the paper is set so that labor is 1/3 in the benchmark case; this normalization is relevant
 *      for the definition of welfare
 * - The technology parameter, which is also left unspecified in the published version is set to in steady state, which 
 *      is the natural normalization also used in the official replication files
 * - For the business cycle moments, the technology shock variance is actually 0.45^2 (not 0.45 as reported in the paper), i.e. 
 *      it is consistent with the number reported in Smets/Wouters (2007). Moreover, business cycle moments rely on a Taylor rule 
 *      with interest rate smoothing of 0, inflation feedback of 1.5, and output feedback of 0.5/4.
 * - As described in footnote 54, the determinacy region in Figure 11 is actually the "determinacy and stability region", i.e. it does not distinguish
 *      whether the Blanchard-Kahn conditions fail because of too many unstable roots (instability, info==3) or because of too few unstable roots (indeterminacy, info==4).
 *      The is replicated in Dynare but not distinguishing between the error codes 3 and 4 returned by resol.m.
 * - Changing the risk aversion parameter from the current log-utility specification requires manually changing the definition of the utility function in 
 *      equation 12 and in the steady_state_model-block. Simply replace the log-utility definition by the commented general CRRA definition.
 *
 * This implementation was written by Johannes Pfeifer. I thank Guido Ascari for providing important clarifications.
 *
 * If you spot mistakes, email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2014-15 Johannes Pfeifer
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


@# define case=1 
%required for TFP shock where three different non-baseline parametrizations are considered

@# define MP_shock=1 
%set to 0 for MP-shock and to 1 for TFP shock

@# define determinacy_plot=1 
%set to 1 if determinacy region should be mapped

var y $y$ //output
    i $i$ //investment
    pi $\pi$ //inflation
    N $N$ //hours worked
    w $w$ //real wage
    p_star ${p^*}$ //target price
    psi $\psi$  //recursive auxiliary variable 1 price setting 
    phi $\phi$ //recursive auxiliary variable 2 price setting 
    A $A$ //TFP
    MC_real $MC$ //real marginal costs 
    real_interest $r$ //real interest rate
    zeta $\zeta$ //preference shock
    s $s$ //price dispersion term
    v $\nu$ //monetary policy shock
    A_tilde ${\tilde A}$ //"effective aggregate productivity"
    Utility $U$ //lifetime utility, recursively defined
    Average_markup 
    Marginal_markup 
    price_adjustment_gap;

varexo e_v e_a e_zeta;

parameters trend_inflation 
    beta $\beta$ //discount factor 
    alpha $\alpha$ //capital share
    phi_par $\varphi$ //Frisch elasticity
    theta $\theta$ //Calvo parameter
    sigma $\sigma$ //Risk aversion
    epsilon $\varepsilon$ //Elasticity of substitution
    Pi_bar ${\bar \pi}$ //gross quarterly steady state inflation
    rho_v ${\rho_\nu}$ //autocorrelation of monetary shock
    rho_a ${\rho_a}$ //autocorrelation of technology shock
    rho_zeta ${\rho_\zeta}$ //autocorrelation of preference shock
    phi_pi ${\phi_\pi}$ //Taylor rule feedback inflation
    phi_y ${\phi_y}$ //Taylor rule output
    Y_bar ${\bar Y}$ //steady state output, set in steady state model block
    var_rho ${\varrho}$ //degree of indexing
    i_bar ${\bar i}$ //steady state interest rate, set in steady state model block
    d_n ${d_n}$ //labor disutility parameter
    rho_i ${\rho_i}$; //interest rate smoothing parameter


%fix labor to 1/3 in zero trend inflation steady state with Frisch elasticity of 1#
beta_ss = 0.99;
alpha_ss = 0;
theta_ss = 0.75;
epsilon_ss = 10;
sigma_ss = 1; %different utility than log case implies that model utility function and its steady state must be manually changed
phi_par_ss=1;
var_rho_ss = 0;
trend_inflation_ss=0;

Pi_bar = (1+0/100)^(trend_inflation_ss/4); %set Pi_bar to reflect quarterly inflation
p_star_ss=((1-theta_ss*Pi_bar^((epsilon_ss-1)*(1-var_rho_ss)))/(1-theta_ss))^(1/(1-epsilon_ss));
s_ss=(1-theta_ss)/(1-theta_ss*Pi_bar^((epsilon_ss*(1-var_rho_ss))/(1-alpha_ss)))*p_star_ss^(-epsilon_ss/(1^-alpha_ss));
N_ss=1/3;
y_ss=(N_ss/s_ss)^(1-alpha_ss);
A_ss=1;
phi_ss=y_ss^(1-sigma_ss)/(1-theta_ss*beta_ss*Pi_bar^((epsilon_ss-1)*(1-var_rho_ss)));
psi_ss=p_star_ss^(1+epsilon_ss*alpha_ss/(1-alpha_ss))*phi_ss/(epsilon_ss/((epsilon_ss-1)*(1-alpha_ss)));
w_ss=psi_ss*(1-theta_ss*beta_ss*Pi_bar^((epsilon_ss*(1-var_rho_ss))/(1-alpha_ss)))/(A_ss^(-1/(1-alpha_ss))*y_ss^(1/(1-alpha_ss)-sigma_ss));
d_n=w_ss/(N_ss^phi_par_ss*y_ss^sigma_ss);


trend_inflation=0; %2, 4, or 6

%set according to FN36
beta = 0.99;
alpha = 0;
theta = 0.75;
epsilon = 10;
sigma = 1; %different utility than log case implies that model utility function and its steady state must be manually changed

rho_v = 0;
rho_a = 0;
rho_zeta = 0;

@# if MP_shock==1
    phi_par = 1; %Frisch elasticity, see FN36
@# else
    @# if case==1
    %analytical case, p. 713
        phi_par = 0; %Frisch elasticity
        rho_a = 0;
    @# endif
    @# if case==2
        phi_par = 0; %Frisch elasticity
        rho_a = 0.95;
    @# endif
    @# if case==3
        phi_par = 3; %Frisch elasticity
        rho_a = 0.95;
    @# endif
@# endif

%Monetary policy according to FN 67
@# if MP_shock==1
phi_pi = 2;
phi_y = 0.5/4;
rho_i=0.8;
@# else
phi_pi = 1.5;
phi_y = 0.5/4;
rho_i=0;
@# endif
        
var_rho = 0;

model;
//1. Euler equation
1/(exp(y)^(sigma)) = beta*(1+exp(i))/(exp(pi(+1))*(exp(y(+1))^(sigma)));
//2. Labor FOC
exp(w) = d_n*exp(zeta)*(exp(N)^phi_par)*(exp(y)^sigma);
//3. Optimal price
exp(p_star) = ((1-theta*(exp(pi(-1))^((1-epsilon)*var_rho))*(exp(pi)^(epsilon-1)))/(1-theta))^(1/(1-epsilon));
//4. FOC price setting
(exp(p_star))^(1+((epsilon*alpha)/(1-alpha))) = (epsilon/((epsilon-1)*(1-alpha)))*exp(psi)/exp(phi);
//5. Recursive LOM price setting for psi
exp(psi) = exp(w)*((exp(A))^(-1/(1-alpha)))*(exp(y)^((1/(1-alpha))-sigma))
            +theta*beta*(exp(pi))^((-var_rho*epsilon)/(1-alpha))*exp(pi(+1))^(epsilon/(1-alpha))*exp(psi(+1));
//6. Recursive LOM price setting for phi
exp(phi) = exp(y)^(1-sigma)+theta*beta*exp(pi)^(var_rho*(1-epsilon))*exp(pi(+1))^(epsilon-1)*exp(phi(+1));
//7. Aggregate production function
exp(N)=exp(s)*(exp(y)/exp(A))^(1/(1-alpha));
//8. LOM price dispersion
exp(s) = (1-theta)*exp(p_star)^(-epsilon/(1-alpha))
        +theta*exp(pi(-1))^((-epsilon*var_rho)/(1-alpha))*exp(pi)^(epsilon/(1-alpha))*exp(s(-1));
//9. Monetary policy rule; reflects FN69
(1+exp(i))/(1+i_bar)=((1+exp(i(-1)))/(1+i_bar))^rho_i*((exp(pi)/Pi_bar)^(phi_pi)*(exp(y)/Y_bar)^(phi_y))^(1-rho_i)*exp(v);

//10. Definition real marginal costs
exp(MC_real)=1/(1-alpha)*exp(w)*exp(A)^(1/(alpha-1))*exp(y)^(alpha/(1-alpha));
//11. Definition real interest rate
exp(real_interest)=(1+exp(i))/(exp(pi(+1)));
//12. Define utility, do not log it as it can be negative; this is the log case
Utility=y-d_n*exp(zeta)*exp(N)^(1+phi_par)/(1+phi_par)+beta*Utility(+1);
// Utility=exp(y)^(1-sigma)/(1-sigma)-d_n*exp(zeta)*exp(N)^(1+phi_par)/(1+phi_par)+beta*Utility(+1);
//13. Monetary shock
v = rho_v*v(-1) + e_v;
//14. Technology shock
A = rho_a*A(-1) + e_a;
//15. Preference shock
zeta = rho_zeta*zeta(-1) + e_zeta;

exp(A_tilde)=exp(A)/exp(s);
exp(Average_markup)=1/exp(MC_real);
exp(Marginal_markup)=exp(p_star)/exp(MC_real);
exp(price_adjustment_gap)=1/exp(p_star);
end;

steady_state_model;
Pi_bar = (1+trend_inflation/100)^(1/4); %set Pi_bar to reflect quarterly inflation
v = 0;
A = 1;
zeta=0;
pi=Pi_bar;
i=1/beta*Pi_bar-1;
i_bar=i;
p_star=((1-theta*Pi_bar^((epsilon-1)*(1-var_rho)))/(1-theta))^(1/(1-epsilon));
s=(1-theta)/(1-theta*Pi_bar^((epsilon*(1-var_rho))/(1-alpha)))*p_star^(-epsilon/(1^-alpha));
y=(p_star^(1+(epsilon*alpha)/(1-alpha))*(epsilon/((epsilon-1)*(1-alpha))*((1-beta*theta*Pi_bar^((epsilon-1)*(1-var_rho)))/(1-beta*theta*Pi_bar^(epsilon*(1-var_rho)/(1-alpha))))*d_n*s^phi_par)^(-1))^((1-alpha)/(phi_par+1));
N=s*y^(1/(1-alpha));
Y_bar=y;
phi=y^(1-sigma)/(1-theta*beta*Pi_bar^((epsilon-1)*(1-var_rho)));
psi=p_star^(1+epsilon*alpha/(1-alpha))*phi/(epsilon/((epsilon-1)*(1-alpha)));
//w=psi*(1-theta*beta*Pi_bar^((epsilon*(1-var_rho))/(1-alpha)))/(A^(-1/(1-alpha))*y^(1/(1-alpha)-sigma));
w=d_n*N^phi_par*y^sigma;
MC_real=1/(1-alpha)*w*A^(1/(alpha-1))*y^(alpha/(1-alpha));
Average_markup=1/MC_real;
Marginal_markup=p_star/MC_real;
real_interest=(1+i)/pi;
price_adjustment_gap=1/p_star;
A_tilde=A/s;
Utility=(1-beta)^(-1)*(log(y)-d_n*N^(1+phi_par)/(1+phi_par));
// Utility=(1-beta)^(-1)*(y^(1-sigma)/(1-sigma)-N^(1-phi_par)/(1+phi_par)); %needed for non-log case
A=log(A);
i=log(i);
p_star=log(p_star);
pi=log(pi);
s=log(s);
y=log(y);
phi=log(phi);
psi=log(psi);
w=log(w);
N=log(N);
MC_real=log(MC_real);
real_interest=log(real_interest);
A_tilde=log(A_tilde);
Average_markup=log(Average_markup);
Marginal_markup=log(Marginal_markup);
price_adjustment_gap=log(price_adjustment_gap);
end;

options_.qz_criterium = 1+1e-6; //make sure the option is set

steady;
check;

@# if MP_shock==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Replicate Figure 14: Impulse Response Functions to a 1 Percent Positive Monetary Policy Shock         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;          
var e_v; stderr 1;   
var e_a; stderr 0;
var e_zeta; stderr 0;
end;

set_param_value('trend_inflation',0)
stoch_simul(irf=7,order=1,irf_plot_threshold=0,nograph,noprint) y pi real_interest s i v; 
irf_0_trend_infl=oo_.irfs;
set_param_value('trend_inflation',2)
stoch_simul(irf=7,order=1,irf_plot_threshold=0,nograph,noprint) y pi real_interest s i v; 
irf_2_trend_infl=oo_.irfs;
set_param_value('trend_inflation',4)
stoch_simul(irf=7,order=1,irf_plot_threshold=0,nograph,noprint) y pi real_interest s i v; 
irf_4_trend_infl=oo_.irfs;
set_param_value('trend_inflation',6)
stoch_simul(irf=7,order=1,irf_plot_threshold=0,nograph,noprint) y pi real_interest s i v; 
irf_6_trend_infl=oo_.irfs;

figure('Name','Figure 14: Impulse Response Functions to a 1 Percent Positive Monetary Policy Shock ')
subplot(2,2,1)
plot(0:options_.irf,[0 irf_0_trend_infl.y_e_v],'k-',0:options_.irf,[0 irf_2_trend_infl.y_e_v],'b--',0:options_.irf,[0 irf_4_trend_infl.y_e_v],'r-.',0:options_.irf,[0 irf_6_trend_infl.y_e_v],'*-')
ylim([-2.5 0.5]); title('Output')
subplot(2,2,2)
plot(0:options_.irf,[0 irf_0_trend_infl.pi_e_v],'k-',0:options_.irf,[0 irf_2_trend_infl.pi_e_v],'b--',0:options_.irf,[0 irf_4_trend_infl.pi_e_v],'r-.',0:options_.irf,[0 irf_6_trend_infl.pi_e_v],'*-')
ylim([-0.80 0]); title('Inflation')
subplot(2,2,3)
plot(0:options_.irf,[0 irf_0_trend_infl.real_interest_e_v],'k-',0:options_.irf,[0 irf_2_trend_infl.real_interest_e_v],'b--',0:options_.irf,[0 irf_4_trend_infl.real_interest_e_v],'r-.',0:options_.irf,[0 irf_6_trend_infl.real_interest_e_v],'*-')
ylim([0 1.4]);title('Real Interest Rate')
subplot(2,2,4)
hh=plot(0:options_.irf,[0 irf_0_trend_infl.s_e_v],'k-',0:options_.irf,[0 irf_2_trend_infl.s_e_v],'b--',0:options_.irf,[0 irf_4_trend_infl.s_e_v],'r-.',0:options_.irf,[0 irf_6_trend_infl.s_e_v],'*-');
ylim([-0.8 0]);title('Price Dispersion')
leg_handle=legend(hh,'pi=0%','pi=2%','pi=4%','pi=6%');

@#else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Replicate Figure 14: Impulse Response Functions to a 1 Percent Positive Technology Shock         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 percent TFP shock
shocks;          
var e_v; stderr 0;   
var e_a; stderr 1;
var e_zeta; stderr 0;
end;

set_param_value('trend_inflation',0)
stoch_simul(irf=14,order=1,nograph,noprint) y pi s ; 
irf_0_trend_infl=oo_.irfs;
set_param_value('trend_inflation',2)
stoch_simul(irf=14,order=1,nograph,noprint) y pi s ; 
irf_2_trend_infl=oo_.irfs;
set_param_value('trend_inflation',4)
stoch_simul(irf=14,order=1,nograph,noprint) y pi s ; 
irf_4_trend_infl=oo_.irfs;

figure('Name','Figure 13: Impulse Response Functions to a 1 Percent Positive Technology Shock')
subplot(3,3,1)
plot(0:options_.irf,[0 irf_0_trend_infl.y_e_a],'k-',0:options_.irf,[0 irf_2_trend_infl.y_e_a],'b--',0:options_.irf,[0 irf_4_trend_infl.y_e_a],'r-.')
ylim([-0.05 0.15]); title('Output')
subplot(3,3,2)
plot(0:options_.irf,[0 irf_0_trend_infl.pi_e_a],'k-',0:options_.irf,[0 irf_2_trend_infl.pi_e_a],'b--',0:options_.irf,[0 irf_4_trend_infl.pi_e_a],'r-.')
ylim([-0.08 0.02]); title('Inflation')
subplot(3,3,3)
plot(0:options_.irf,[0 irf_0_trend_infl.s_e_a],'k-',0:options_.irf,[0 irf_2_trend_infl.s_e_a],'b--',0:options_.irf,[0 irf_4_trend_infl.s_e_a],'r-.')
ylim([-0.02 0.005]);title('Price Dispersion')


%Consider case 2 with rho_a=0.95
set_param_value('rho_a',0.95)
set_param_value('trend_inflation',0)
stoch_simul(irf=14,order=1,nograph,noprint) y pi s ; 
irf_0_trend_infl=oo_.irfs;
set_param_value('trend_inflation',2)
stoch_simul(irf=14,order=1,nograph,noprint) y pi s ; 
irf_2_trend_infl=oo_.irfs;
set_param_value('trend_inflation',4)
stoch_simul(irf=14,order=1,nograph,noprint) y pi s ; 
irf_4_trend_infl=oo_.irfs;
subplot(3,3,4)
plot(0:options_.irf,[0 irf_0_trend_infl.y_e_a],'k-',0:options_.irf,[0 irf_2_trend_infl.y_e_a],'b--',0:options_.irf,[0 irf_4_trend_infl.y_e_a],'r-.')
ylim([0 1.5])
subplot(3,3,5)
plot(0:options_.irf,[0 irf_0_trend_infl.pi_e_a],'k-',0:options_.irf,[0 irf_2_trend_infl.pi_e_a],'b--',0:options_.irf,[0 irf_4_trend_infl.pi_e_a],'r-.')
ylim([-0.4 0])
subplot(3,3,6)
plot(0:options_.irf,[0 irf_0_trend_infl.s_e_a],'k-',0:options_.irf,[0 irf_2_trend_infl.s_e_a],'b--',0:options_.irf,[0 irf_4_trend_infl.s_e_a],'r-.')
ylim([-0.65 0.2])

%Consider case 3 with Frisch elasticity of 3
set_param_value('phi_par',3)
set_param_value('trend_inflation',0)
stoch_simul(irf=14,order=1,nograph,noprint) y pi s ; 
irf_0_trend_infl=oo_.irfs;
set_param_value('trend_inflation',2)
stoch_simul(irf=14,order=1,nograph,noprint) y pi s ; 
irf_2_trend_infl=oo_.irfs;
set_param_value('trend_inflation',4)
stoch_simul(irf=14,order=1,nograph,noprint) y pi s ; 
irf_4_trend_infl=oo_.irfs;

subplot(3,3,7)
plot(0:options_.irf,[0 irf_0_trend_infl.y_e_a],'k-',0:options_.irf,[0 irf_2_trend_infl.y_e_a],'b--',0:options_.irf,[0 irf_4_trend_infl.y_e_a],'r-.')
ylim([0 1.55])
subplot(3,3,8)
plot(0:options_.irf,[0 irf_0_trend_infl.pi_e_a],'k-',0:options_.irf,[0 irf_2_trend_infl.pi_e_a],'b--',0:options_.irf,[0 irf_4_trend_infl.pi_e_a],'r-.')
ylim([-0.4 0])
subplot(3,3,9)
hh=plot(0:options_.irf,[0 irf_0_trend_infl.s_e_a],'k-',0:options_.irf,[0 irf_2_trend_infl.s_e_a],'b--',0:options_.irf,[0 irf_4_trend_infl.s_e_a],'r-.');
ylim([-1 0.5])
leg_handle=legend(hh,'pi=0%','pi=2%','pi=4%');
@#endif

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Replicate business cycle moments reported on page 717  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shocks;          
var e_v; stderr 0;   
var e_a=0.45^2;
var e_zeta; stderr 0;
end;

set_param_value('trend_inflation',0);
set_param_value('phi_par',1);
set_param_value('rho_a',0.95);
set_param_value('rho_i',0); % was set to 0 for moments reported
set_param_value('phi_pi',1.5);
set_param_value('phi_y',0.5/4);
stoch_simul(irf=0,order=1,nograph,noprint) y pi; 
output_var_low_target=oo_.var(strmatch('y',var_list_,'exact'),strmatch('y',var_list_,'exact'));
inflation_var_low_target=oo_.var(strmatch('pi',var_list_,'exact'),strmatch('pi',var_list_,'exact'));
set_param_value('trend_inflation',4);
stoch_simul(irf=0,order=1,nograph,noprint) y pi; 
output_var_high_target=oo_.var(strmatch('y',var_list_,'exact'),strmatch('y',var_list_,'exact'));
inflation_var_high_target=oo_.var(strmatch('pi',var_list_,'exact'),strmatch('pi',var_list_,'exact'));
fprintf('Output Standard Deviation: \t %4.3f \t %4.3f\n',sqrt(output_var_low_target),sqrt(output_var_high_target));
fprintf('Inflation Standard Deviation: \t %4.3f \t %4.3f\n',sqrt(inflation_var_low_target),sqrt(inflation_var_high_target))

verbatim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Generate Figure 7: The Cost of Price Dispersion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Name','The Cost of Price Dispersion')
    trend_inflation_vector=[2,4];
    theta_vector=linspace(0.5,0.87,50);
    output_loss=NaN(length(theta_vector),2);
    for trend_inflation_iter=1:2
        set_param_value('trend_inflation',trend_inflation_vector(trend_inflation_iter));
        for iter=1:length(theta_vector)
            set_param_value('theta',theta_vector(iter));
            steady;
            output_loss(iter,trend_inflation_iter)=oo_.steady_state(strmatch('A_tilde',M_.endo_names,'exact'))*100;
        end
    end    
    subplot(1,3,1)
    plot(theta_vector,output_loss(:,1)-output_loss(1,1),'-',theta_vector,output_loss(:,2)-output_loss(1,2),'--'); %normalize relative to first value
    xlabel('theta')
    
    set_param_value('theta',0.75);
    trend_inflation_vector=[2,4];
    epsilon_vector=linspace(1.01,14,50);
    output_loss=NaN(length(epsilon_vector),2);
    for trend_inflation_iter=1:2
        set_param_value('trend_inflation',trend_inflation_vector(trend_inflation_iter));
        for iter=1:length(epsilon_vector)
            set_param_value('epsilon',epsilon_vector(iter));
            steady;
            output_loss(iter,trend_inflation_iter)=oo_.steady_state(strmatch('A_tilde',M_.endo_names,'exact'))*100;%(1/exp(oo_.steady_state(strmatch('s',M_.endo_names,'exact')))-1)*100;
        end
    end    
    subplot(1,3,2)
    plot(epsilon_vector,output_loss(:,1)-output_loss(1,1),'-',epsilon_vector,output_loss(:,2)-output_loss(1,2),'--'); %normalize relative to first value
    xlabel('epsilon')
    
    set_param_value('theta',0.75); %reset to baseline
    set_param_value('epsilon',10); %reset to baseline
    trend_inflation_vector=0:0.5:8;
    output_loss=NaN(length(trend_inflation_vector),1);
    for iter=1:length(trend_inflation_vector)
        set_param_value('trend_inflation',trend_inflation_vector(iter));
        steady;
        output_loss(iter,1)=oo_.steady_state(strmatch('A_tilde',M_.endo_names,'exact'))*100;
    end
    subplot(1,3,3)
    plot(trend_inflation_vector,output_loss)
    xlabel('Trend Inflation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Generate Figure 8: Trend Inflation and Steady State Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    figure('Name','Trend Inflation and Steady State Variables')
    trend_inflation_vector=0:0.5:8;
    utility=NaN(length(trend_inflation_vector),1);
    output=NaN(length(trend_inflation_vector),1);
    marg_markup=NaN(length(trend_inflation_vector),1);
    ave_markup=NaN(length(trend_inflation_vector),1);
    price_adjust_gap=NaN(length(trend_inflation_vector),1);
    for iter=1:length(trend_inflation_vector)
        set_param_value('trend_inflation',trend_inflation_vector(iter));
        steady;
        utility(iter,1)=oo_.steady_state(strmatch('Utility',M_.endo_names,'exact'));
        output(iter,1)=oo_.steady_state(strmatch('y',M_.endo_names,'exact'));
        marg_markup(iter,1)=oo_.steady_state(strmatch('Marginal_markup',M_.endo_names,'exact'));
        ave_markup(iter,1)=oo_.steady_state(strmatch('Average_markup',M_.endo_names,'exact'));
        marg_markup_eq_37(iter,1)=log(epsilon/(epsilon-1)*(1-beta*theta*(M_.params(strmatch('Pi_bar',M_.param_names,'exact')))^(epsilon-1))/(1-beta*theta*(M_.params(strmatch('Pi_bar',M_.param_names,'exact')))^(epsilon)));
        price_adjust_gap(iter,1)=oo_.steady_state(strmatch('price_adjustment_gap',M_.endo_names,'exact'));
    end
    if max(abs(marg_markup-marg_markup_eq_37))>1e-8
        error('Wrong results')
    end
    subplot(1,3,1)
    plot(trend_inflation_vector,(output-output(1,1))*100)
    xlabel('Trend Inflation')
    ylabel('Steady state output')
    subplot(1,3,2)
    plot(trend_inflation_vector,(ave_markup-ave_markup(1,1))*100,'-',trend_inflation_vector,(marg_markup-marg_markup(1,1))*100,'--',trend_inflation_vector,(price_adjust_gap-price_adjust_gap(1,1))*100,'.')
    xlabel('Trend Inflation')
    ylabel('SS price adj. gap, marg. and ave. markup')
    subplot(1,3,3)
    plot(trend_inflation_vector,(utility-utility(1,1))./abs(utility(1,1))*100)
    xlabel('Trend Inflation')
    ylabel('Steady state welfare')
    
end;
    
set_param_value('phi_par',1); %reset to baseline
set_param_value('rho_i',0); %reset to baseline
set_param_value('theta',0.75); %reset to baseline
set_param_value('epsilon',10); %reset to baseline

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Generate Figure 11: The Determinacy Region and Trend Inflation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@# if determinacy_plot==1
    phi_pi_vec=linspace(0,5,200);
    phi_y_vec=linspace(-1,5,200);
    [phi_pi_mat,phi_y_mat]=meshgrid(phi_pi_vec,phi_y_vec);
    trend_inflation_vector=[0,2,4,6,8];
    Z_plot_total=zeros(size(phi_pi_mat));
    for trend_inflation_iter=1:length(trend_inflation_vector)
        set_param_value('trend_inflation',trend_inflation_vector(trend_inflation_iter));
        info_mat=NaN(size(phi_pi_mat));
        for phi_pi_iter=1:length(phi_pi_vec)
            for phi_y_iter=1:length(phi_y_vec)
                set_param_value('phi_pi',phi_pi_mat(phi_pi_iter,phi_y_iter));
                set_param_value('phi_y',phi_y_mat(phi_pi_iter,phi_y_iter));
                [dr,info]=resol(0,M_,options_,oo_);
                info_mat(phi_pi_iter,phi_y_iter)=info(1);
            end
        end
        Z_plot=zeros(size(info_mat));
        Z_plot(info_mat==0)=1;
        figure
        contourf(phi_pi_mat,phi_y_mat,Z_plot,1)
        Z_plot_total(info_mat==0)=trend_inflation_iter;
        xlabel('\phi_\pi')
        ylabel('\phi_y')
        title([num2str(trend_inflation_vector(trend_inflation_iter)) '%'])
    end    

    figure('Name','Figure 11: The Determinacy Region and Trend Inflation')
    contourf(phi_pi_mat,phi_y_mat,Z_plot_total,5);
    colormap(hot);
@#endif
                
// A similar figure could have been obtained using Dynare's dynare_sensitivity command. The difference
// is that the parameter would have been randomly sampled from the prior instead of uniformly on the grid
// Moreover, many draws (nsam) are required to get a mapping figure of the type obtained above.
//         
// %specify parameters for which to map sensitivity
// estimated_params;
// phi_pi,uniform_pdf,(0+6)/2,sqrt(12)^(-1)*(6-0);
// phi_y,uniform_pdf,((-1)+6)/2,sqrt(12)^(-1)*(6-(-1));      
// end;
// 
// varobs y pi; //some observables must be specified for sensitivity command, inessential for results
// options_.nograph=0; %enable graphs again
// 
// set_param_value('trend_inflation',2);
// dynare_sensitivity(prior_range=0,stab=1,nsam=5000);
