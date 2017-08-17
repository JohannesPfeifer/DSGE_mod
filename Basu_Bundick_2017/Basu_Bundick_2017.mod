/* Replicates the IRFs at the stochastic steady/ergodic mean in the absence of shocks 
 * by Basu/Bundick (2017): "Uncertainty shocks in a model of effective demand", 
 * Econometrica, 85(3), pp. 937-958
 *
 * Notes:
 *  - Due to pruning, one cannot simply use the stochastic steady state as the starting point 
 *      for simult_.m as one would do with the deterministic steady state in a linear model. The reason
 *      is the pruned state space where one would need to provide the first, second, and third order components
 *      of the stochastic steady state. Dynare does currently not yet support this. For this reason, the impulse
 *      period is simply appended to the simulation for computing the stochastic steady state
 *  - Basu/Bundick only use 200 periods of burn-in for computing the stochastic steady state. This is not sufficient
 *      for convergence as can be easily verified by plotting a longer simulation. Setting true_stochastic_steady_state_IRFs=1
 *      therefore uses a longer burn-in, leading to slightly different IRFs.
*/
/*
 * Copyright (C) 2016-17 Benjamin Born and Johannes Pfeifer
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

@#define sticky_prices = 1
@#define CRRA = 0
@#define true_stochastic_steady_state_IRFs = 1

var C       ${C}$                               (long_name='consumption')
    N       ${N}$                               (long_name='labor')
    P_E     ${\left({\frac{P^E}{P}}\right)}$ 	(long_name='Real price of shares')
    D_E     ${\left({\frac{D^E}{P}}\right)}$ 	(long_name='Real Dividends of shares')
    R_R     ${R^R}$                             (long_name='risk-free rate')
    W       ${\left({\frac{W}{P}}\right)}$      (long_name='real wage')
    V       ${V}$                               (long_name='Value function')
    E_t_V_tp1_1_minus_sigma ${\left({E_t V_{t+1}^{1-\sigma}}\right)}$   (long_name='Auxiliary Variable for EZ')
    Y       $Y$                                 (long_name='Output')
    K       $K$                                 (long_name='capital')
    I       $I$                                 (long_name='Investment')
    D       ${\frac{D}{P}}$                     (long_name='Cash flows')
    Pi      $\Pi$                               (long_name='inflation')
    R_K     ${R^k}$                             (long_name='marginal revenue product of capital')
    q       $q$                                 (long_name='Tobins q')
    mu      $\mu$                               (long_name='markup')
    sigma_a ${\sigma_a}$                        (long_name='preference shock volatility')
    a       $a$                                 (long_name='preference shock level')
    Z       $Z$                                 (long_name='technology shock level')
    R       $R$                                 (long_name='policy rate')
    Xi      $\Xi$                               (long_name='marginal costs')
    u       $u$                                 (long_name='capacity utilization')
    PHI     ${\Phi}$                            (long_name='Rotemberg price adjustment costs')
    R_E     ${R^E}$                             (long_name='return to equity')
    E_R_E   ${E_t(R^E_{t+1})}$                  (long_name='expected return to equity')
    E_R_E_squared   ${(E_t(R^E_{t+1})^2)}$      (long_name='expected squared return to equity')
    cond_var_R_E    ${VAR_t(R^E_{t+1})}$        (long_name='conditional variance of return to equity')
    E_M_tp1         ${(E_t(M_{t+1}))}$          (long_name='expected SDF')
    E_R_E_risk_neutral ${(E_t(R^E_{t+1})^2)}$   (long_name='expected SDF squared')
    E_M_tp1_R_E     ${E_t(M_{t+1} R^E_{t+1})}$  (long_name='expecation of SDF times equity return')
    E_M_tp1_R_E_squared ${E_t(M_{t+1} (R^E_{t+1})^2)}$  (long_name='expecation of SDF times squared equity return')
    cond_var_R_E_risk_neutral   ${VAR_t(R^E_{t+1})}$    (long_name='conditional variance of risk-neutral return to equity')
    log_Y               ${\hat Y}$              (long_name='log output')
    log_C               ${\hat C}$              (long_name='log consumption')
    log_I               ${\hat I}$              (long_name='log investment')
    log_mu              ${\hat \mu}$            (long_name='log markup')
    log_N               ${\hat N}$              (long_name='log hours')
    log_sigma_a         ${\hat {\sigma_a}}$     (long_name='log preference shock volatility')
    log_W               ${\hat W}$              (long_name='log wage')
    pi_annualized       ${\Pi^{ann}}$           (long_name='annualized inflation')
    R_annualized        ${\R^{ann}}$            (long_name='annualized interest rate')
    R_K_annualized      ${\R^{K,ann}}$          (long_name='annualized capital return')
    R_R_annualized      ${\R^{R,ann}}$          (long_name='annualized risk-free rate')
    log_pi_annualized   ${\log(\Pi^{ann})}$     (long_name='annualized log inflation')
    log_R_annualized    ${\log(\R^{ann})}$      (long_name='annualized log interest rate')
    log_R_K_annualized  ${\log(\R^{K,ann})}$    (long_name='annualized log capital return')
    log_R_R_annualized  ${\log(\R^{R,ann})}$    (long_name='annualized log risk-free interest rate')
;

predetermined_variables K;

varexo eps_a        ${\varepsilon^a}$           (long_name='preference shock')
    eps_sigma_a     ${\varepsilon^{\sigma^a}}$  (long_name='preference shock volatility')
    eps_z           ${\varepsilon^z}$       	(long_name='technology shock')
    eps_m           ${\varepsilon^m}$       	(long_name='monetary policy shock')
;

parameters siggma   ${\sigma}$          (long_name='risk aversion')
    theta_v         ${\theta_v}$        (long_name='uncertainty resolution preference')
    psii            $\psi$              (long_name='intertemporal elasticity of substitution')
    eta             $\eta$              (long_name='share of consumption in Cobb-Douglas aggregator')
    betta           $\beta$             (long_name='discount factor')
    theta_mu        ${\theta_\mu}$      (long_name='demand elasticity')
    alppha          $\alpha$            (long_name='labor share')
    Psi             $\Psi$              (long_name='fixed costs')
    delta_0         $\delta_0$          (long_name='ss depreciation')
    delta_1         $\delta_1$          (long_name='variable cu')
    delta_2         $\delta_2$          (long_name='variable cu')
    phi_k           $\phi_k$            (long_name='capital adjustment costs')
    phi_p           $\phi_P$            (long_name='price adjustment costs')
    Pi_bar          ${\bar \Pi}$        (long_name='steady state inflation')
    nu              $\nu$               (long_name='share of bonds in capital')
    rho_r           $\rho_r$            (long_name='interest smoothing')
    rho_pi          $\rho_\pi$          (long_name='inflation feedback')
    rho_y           $\rho_y$            (long_name='output growth feedback')
    log_R_bar       $\bar r$            (long_name='SS net interest rate')
    sigma_a_bar     ${\sigma^a}$        (long_name='SS volatility preference shock')
    sigma_sigma_a   ${\sigma^{\sigma^a}}$     (long_name='SS volatility of uncertainty shock')
    rho_a           $\rho_a$            (long_name='persistence preference shock')
    rho_sigma_a     $\rho_{\sigma^a}$   (long_name='pseristence volatility shock')
    a_bar           $a$                 (long_name='SS preference shock')
    sigma_z_bar     ${\sigma^z}$        (long_name='SS volatility technology shock')
    rho_z           $\rho_z$            (long_name='persistence technology shock')
    PF_normalization $A$                (long_name='Normalization production function')
    V_normalization ${V^{norm}}$        (long_name='Normalization value function')
    Frisch_target   ${Frisch elast.}$   (long_name='Target Frisch elasticity')
    sigma_eps_m     ${\sigma_m}$        (long_name='monetary policy shock volatility')
;


alppha =0.333;
betta =0.994;
delta_0 =0.025;
delta_1 =0; %set in SS file 
delta_2 =0.00031;
@#if sticky_prices == 1
    phi_p = 100;
@#else
    phi_p =0;
@# endif

Pi_bar =1.005;
log_R_bar=log(Pi_bar/betta);
rho_r =0;
rho_pi =1.5;
rho_y =0.2;
@#if CRRA == 1
siggma = 1.1;
psii = 1/siggma;
@#else
siggma = 80;
psii = 0.95;
@#endif
theta_mu = 6;
a_bar =1;
Frisch_target=2;

phi_k =2.0901;

rho_a =0.93564;
sigma_a_bar =0.0026251;
rho_sigma_a =0.74227;
sigma_sigma_a =0.0025022;

rho_z =0.98793;
sigma_z_bar =0.0012857;
nu=0.9;     %irrelevant due to Modigliani-Miller
sigma_eps_m=0;

eta=0;      %set in SS file 
Psi = 0;    %fixed costs, set in SS file

model;
@#if CRRA == 1
    [name='Definition value function']
    V=V_normalization*a*(C^eta*(1-N)^(1-eta))^(1-siggma)+betta*V(+1);
    #M=betta*(a(+1)/a)*((C(+1)^eta*(1-N(+1))^(1-eta))/(C^eta*(1-N)^(1-eta)))^(1-siggma)*(C/C(+1));
@#else
    [name='Definition value function']
    V=(V_normalization*a*(C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v) +betta*(E_t_V_tp1_1_minus_sigma)^(1/theta_v))^(theta_v/(1-siggma));
    #M=betta*(a(+1)/a)*((C(+1)^eta*(1-N(+1))^(1-eta))/(C^eta*(1-N)^(1-eta)))^((1-siggma)/theta_v)*(C/C(+1))*(V(+1)^(1-siggma)/E_t_V_tp1_1_minus_sigma)^(1-1/theta_v);
@#endif
[name='Auxiliary variable needed for Epstein-Zin preferences']
E_t_V_tp1_1_minus_sigma=V(+1)^(1-siggma);
[name='Budget constraint']
C+P_E+1/R_R*nu*K(+1)=W*N+(D_E+P_E)+nu*K;
[name='(10) labor FOC']
(1-eta)/eta*C/(1-N)=W;
[name='(11) Stock FOC']
P_E=M*(D_E(+1)+P_E(+1));
[name='(12) Firm Bond FOC']
1=R_R*M;
[name='Production function']
Y=PF_normalization*(u*K)^alppha*(Z*N)^(1-alppha)-Psi;
[name='LOM Capital']
K(+1)=(1-(delta_0 + delta_1*(u-1)+delta_2/2*(u-1)^2)-phi_k/2*(I/K-delta_0)^2)*K+I;
[name='Cash flows']
D=Y-W*N-I-phi_p/2*(Pi/Pi_bar-1)^2*Y;
[name='Labor FOC']
W*N=(1-alppha)*Xi*PF_normalization*(u*K)^alppha*(Z*N)^(1-alppha);
[name='Capital FOC']
R_K*(u*K)=alppha*Xi*PF_normalization*(u*K)^alppha*(Z*N)^(1-alppha);
[name='Utilization FOC']
q*(delta_1+delta_2*(u-1))*u*K=alppha*Xi*PF_normalization*(u*K)^alppha*(Z*N)^(1-alppha);
[name='Pricing FOC']
phi_p*(Pi/Pi_bar-1)*(Pi/Pi_bar)=(1-theta_mu)+theta_mu*Xi+phi_p*
    M*(Y(+1)/Y)*(Pi(+1)/Pi_bar-1)*(Pi(+1)/Pi_bar);
[name='FOC capital']
q=M*(u(+1)*R_K(+1)+q(+1)*(1-(delta_0 + delta_1*(u(+1)-1)+delta_2/2*(u(+1)-1)^2)-phi_k/2*(I(+1)/K(+1)-delta_0)^2 +phi_k*(I(+1)/K(+1)-delta_0)*(I(+1)/K(+1))));
[name='FOC investment']
1/q=1-phi_k*(I/K-delta_0);
[name='Dividends']
D_E=D-nu*(K-1/R_R*K(+1));
[name='Taylor Rule']
log(R)=rho_r*log(R(-1))+(1-rho_r)*(log_R_bar+rho_pi*(log(Pi)-log(Pi_bar))+rho_y*(log(Y)-log(Y(-1))))+sigma_eps_m*eps_m;
[name='Bond FOC']
1=betta*(a(+1)/a)*((C(+1)^eta*(1-N(+1))^(1-eta))/(C^eta*(1-N)^(1-eta)))^((1-siggma)/theta_v)*(C/C(+1))*(V(+1)^(1-siggma)/E_t_V_tp1_1_minus_sigma)^(1-1/theta_v)*(R/Pi(+1));
[name='Markup definition']
mu=1/Xi;
[name='Rotemberg cost definition']
PHI=1+phi_p/2*(Pi/Pi_bar-1)^2*Y;
[name='Definition equity return']
R_E=(D_E + P_E)/P_E(-1);
[name='Definition expected equity return']
E_R_E=R_E(+1);
[name='Definition expected squared equity return']
E_R_E_squared=R_E(+1)^2;
[name='Definition conditional variance equity return']
cond_var_R_E=E_R_E_squared-E_R_E^2;

E_M_tp1=M;
E_R_E_risk_neutral=R_E(+1)/E_M_tp1;
E_M_tp1_R_E=M*R_E(+1);
E_M_tp1_R_E_squared=M*R_E(+1)^2;
cond_var_R_E_risk_neutral=E_M_tp1_R_E_squared/E_M_tp1-(E_M_tp1_R_E/E_M_tp1)^2;


[name='Preference shock process level equation']
a=(1-rho_a)*a_bar+rho_a*a(-1)+sigma_a(-1)*eps_a;
[name='Preference shock volatility process']
sigma_a=(1-rho_sigma_a)*sigma_a_bar+rho_sigma_a*sigma_a(-1)+sigma_sigma_a*eps_sigma_a;
[name='TFP shock process']
Z=(1-rho_z)*1+rho_z*Z(-1)-sigma_z_bar*eps_z;

%Observed variables
[name='Definition log output']
log_Y=log(Y);
[name='Definition log consumption']
log_C=log(C);
[name='Definition log investment']
log_I=log(I);
[name='Definition log markup']
log_mu=log(mu);
[name='Definition log hours']
log_N=log(N);
[name='Definition log volatility']
log_sigma_a=log(sigma_a);
[name='Definition log wage']
log_W =log(W);
[name='Definition annualized inflation']
pi_annualized=1+(4*(Pi-1));
[name='Definition annualized interest rate']
R_annualized=1+(4*(R-1));
[name='Definition annualized return to capital']
R_K_annualized=4*R_K;
[name='Definition annualized risk-free rate']
R_R_annualized=1+(4*(R_R-1));
[name='Definition annualized log inflation']
log_pi_annualized=log(pi_annualized);
[name='Definition annualized log interest rate']
log_R_annualized=log(R_annualized);
[name='Definition annualized log return to capital']
log_R_K_annualized=log(R_K_annualized);
[name='Definition annualized log risk-free rate']
log_R_R_annualized=log(R_R_annualized);
end;


write_latex_dynamic_model;
write_latex_static_model;
write_latex_definitions;

shocks;
var eps_a; stderr 1;
var eps_sigma_a; stderr 1;
var eps_z; stderr 1;
var eps_m; stderr 0;
end;

steady;

options_.TeX=1;
write_latex_dynamic_model;
write_latex_parameter_table;

stoch_simul(order=3,pruning,k_order_solver,noprint,irf=0
            ) Y log_Y log_C log_I log_mu log_N log_sigma_a log_pi_annualized log_R_annualized log_R_R_annualized log_W log_R_K_annualized cond_var_R_E Z a sigma_a;

y_pos 	= strmatch('Y',M_.endo_names,'exact');
c_pos 	= strmatch('C',M_.endo_names,'exact');
inv_pos = strmatch('I',M_.endo_names,'exact');
mu_pos 	= strmatch('mu',M_.endo_names,'exact');
n_pos 	= strmatch('N',M_.endo_names,'exact');
w_pos 	= strmatch('W',M_.endo_names,'exact');
pi_pos 	= strmatch('Pi',M_.endo_names,'exact');
r_pos 	= strmatch('R',M_.endo_names,'exact');
R_R_pos = strmatch('R_R',M_.endo_names,'exact');
sigma_a_pos = strmatch('sigma_a',M_.endo_names,'exact');
Z_pos    = strmatch('Z',M_.endo_names,'exact');
conditional_variance_R_E_pos = strmatch('cond_var_R_E',M_.endo_names,'exact');


IRF_periods=20;
@#if true_stochastic_steady_state_IRFs
    burnin=5000; %periods for convergence
@#else
    burnin=5000; %periods for convergence
@#endif

shock_mat_with_zeros=zeros(burnin+IRF_periods,M_.exo_nbr); %shocks set to 0 to simulate without uncertainty
IRF_no_shock_mat = simult_(oo_.dr.ys,oo_.dr,shock_mat_with_zeros,options_.order)'; %simulate series
stochastic_steady_state=IRF_no_shock_mat(1+burnin,:); % stochastic_steady_state/EMAS is any of the final points after burnin

shock_mat = zeros(burnin+IRF_periods,M_.exo_nbr);
shock_mat(1+burnin,strmatch('eps_sigma_a',M_.exo_names,'exact'))= 1;
IRF_mat = simult_(oo_.dr.ys,oo_.dr,shock_mat,options_.order)';

IRF_mat_percent_from_SSS = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,:)-IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,:))./repmat(stochastic_steady_state,IRF_periods,1); %only valid for variables not yet logged

%scale IRFs as reqired
y_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,y_pos);
c_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,c_pos);
inv_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,inv_pos);
mu_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,mu_pos);
n_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,n_pos);
w_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,w_pos);
pi_vola_IRF 	= 400*IRF_mat_percent_from_SSS(:,pi_pos);
r_vola_IRF 		= 400*IRF_mat_percent_from_SSS(:,r_pos);
R_R_vola_IRF 	= 400*IRF_mat_percent_from_SSS(:,R_R_pos);
sigma_a_vola_IRF= 100*IRF_mat_percent_from_SSS(:,sigma_a_pos);
Z_vola_IRF      = 100*IRF_mat_percent_from_SSS(:,Z_pos);

vxo_shock 		= 100*sqrt(4*(IRF_mat(burnin+1:burnin+IRF_periods,conditional_variance_R_E_pos)));
vxo_mean=100*sqrt(4*(stochastic_steady_state(conditional_variance_R_E_pos)));
vxo_vola_IRF 	= 100*(vxo_shock/vxo_mean-1);


hh=figure;
figure(hh)   
subplot(3,3,1)
hold on
plot(y_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Output','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,2)
hold on
plot(c_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Consumption','FontSize',14)
ylabel('Percent','FontSize',12)
%ylim([-0.3 0.1]);set(gca,'YTick',[-0.3:0.1:0.1],'FontSize',12);


figure(hh)   
subplot(3,3,3)
hold on
plot(inv_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Investment','FontSize',14)
ylabel('Percent','FontSize',12)
%ylim([-0.6 0.4]);set(gca,'YTick',[-0.6:0.2:0.4],'FontSize',12);

figure(hh)   
subplot(3,3,4)
hold on
plot(mu_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Markup','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,5)
hold on
plot(n_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Hours Worked','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,6)
hold on
plot(w_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Real Wage','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,7)
hold on
plot(R_R_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Real Interest Rate','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,8)
hold on
plot(pi_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Inflation','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,9)
hold on
plot(sigma_a_vola_IRF,'b-','LineWidth',3)
plot(sigma_a_vola_IRF,'r--','LineWidth',3,'HandleVisibility','off')
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Preference Shock Volatility','FontSize',14)
ylabel('Percent','FontSize',12)

figure
hold on
plot(vxo_shock,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('VXO','FontSize',14)
ylabel('Percent','FontSize',12)