/*
 * This file implements the baseline New Keynesian model of Jordi Galí (2015): Monetary Policy, Inflation,
 * and the Business Cycle, Princeton University Press, Second Edition, Chapter 3
 *
 * Note that this mod-file implements the non-linear first order conditions and that the IRFs show the log-deviations
 * from steady state.
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER
 *
 * Notes:
 *  - in the LOM for the discount rate shock z the shock enters with a minus sign in this mod-file to generate the 
 *      IRF to a -0.5% shock
 *  - the IRF for the nominal rate in Figure 3.6 "Dynamic Responses to a Technology Shock: Money Supply Rule", p. 81
 *      is wrong. It should be identically 0 as can be seen in this mod-file and Galí's slide set accompanying this chapter
 *  - For nonlinear models, the proper way to specify standard deviations is in the form 0.01 for 1%. One must not multiply by 
 *      100 here, because the uncertainty correct will be wrong in this case.
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model. 
 */

/*
 * Copyright (C) 2016 Johannes Pfeifer
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

%define whether to use interest rate or money growth rate rule 
@#define money_growth_rule=1

var C               ${C}$           (long_name='Consumption')
    W_real          ${\frac{W}{P}}$ (long_name='Real Wage')
    Pi              ${\Pi}$         (long_name='inflation')
    A               ${A}$           (long_name='AR(1) technology process')
    N               ${N}$           (long_name='Hours worked')
    R               ${R^n}$         (long_name='Nominal Interest Rate') 
    realinterest    ${R^{r}}$       (long_name='Real Interest Rate')
    Y               ${Y}$           (long_name='Output') 
    Q               ${Q}$           (long_name='Bond price')
    Z               ${Z}$           (long_name='AR(1) preference shock process')
    S               ${S}$           (long_name='Price dispersion')
    Pi_star         ${\Pi^*}$       (long_name='Optimal reset price')
    x_aux_1         ${x_1}$         (long_name='aux. var. 1 recursive price setting')
    x_aux_2         ${x_2}$         (long_name='aux. var. 2 recursive price setting')
    MC              ${mc}$          (long_name='real marginal costs')
    M_real          ${M/P}$         (long_name='real money stock')
    i_ann           ${i^{ann}}$     (long_name='annualized nominal interest rate')
    pi_ann          ${\pi^{ann}}$   (long_name='annualized inflation rate')
    r_real_ann      ${r^{r,ann}}$   (long_name='annualized real interest rate')
    P               ${P}$           (long_name='price level')
    log_m_nominal   ${log(M)}$      (long_name='log nominal money stock')
    log_y           ${log(Y)}$      (long_name='log output')
    log_W_real      ${log(W/P)}$    (long_name='log real wage')
    log_N           ${log(N)}$      (long_name='log hours')
    log_P           ${log(P)}$      (long_name='log price level')
    log_A           ${log(A)}$      (long_name='log technology level')
    log_Z           ${log(Z)}$      (long_name='log preference shock')
    @#if money_growth_rule==0
        nu                  ${\nu}$             (long_name='AR(1) monetary policy shock process')    
    @#else
        money_growth        ${\Delta m_q}$      (long_name='money growth')
        money_growth_ann    ${\Delta m^{ann}}$  (long_name='money growth annualized')
    @#endif    
    ;

varexo eps_a        ${\varepsilon_a}$   (long_name='technology shock')
       eps_z        ${\varepsilon_z}$   (long_name='preference shock')
       @#if money_growth_rule==0
            eps_nu  ${\varepsilon_\nu}$ (long_name='monetary policy shock')
       @#else   
           eps_m    ${\varepsilon_m}$   (long_name='money supply shock innovation')
       @#endif
       ;   

parameters alppha       ${\alpha}$      (long_name='capital share')
    betta               ${\beta}$       (long_name='discount factor')
    rho_a               ${\rho_a}$      (long_name='autocorrelation technology shock')
    @#if money_growth_rule==0
        rho_nu          ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')
    @#else   
        rho_m           ${\rho_{\zeta}}$ (long_name='autocorrelation monetary demand shock')
    @#endif
    rho_z               ${\rho_{z}}$    (long_name='autocorrelation monetary demand shock')
    siggma              ${\sigma}$      (long_name='inverse EIS')
    varphi              ${\varphi}$     (long_name='inverse Frisch elasticity')
    phi_pi              ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    phi_y               ${\phi_{y}}$    (long_name='output feedback Taylor Rule')
    eta                 ${\eta}$        (long_name='semi-elasticity of money demand')
    epsilon             ${\epsilon}$    (long_name='demand elasticity')
    theta               ${\theta}$      (long_name='Calvo parameter')
    tau                 ${\tau}$      (long_name='labor subsidy')
    ;
    
%----------------------------------------------------------------
% Parametrization, p. 67  and p. 113-115
%----------------------------------------------------------------
siggma = 1;
varphi=5;
phi_pi = 1.5;
phi_y  = 0.125;
theta=3/4;
@#if money_growth_rule==0
    rho_nu =0.5;
@#else   
    rho_m=0.5; %footnote 11, p. 115
@#endif
rho_z  = 0.5;
rho_a  = 0.9;
betta  = 0.99;
eta  =3.77; %footnote 11, p. 115
alppha=1/4;
epsilon=9;
tau=0; //1/epsilon;

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model;
    [name='FOC Wages, eq. (2)']
    W_real=C^siggma*N^varphi;
    [name='Euler equation eq. (3)']
    Q=betta*(C(+1)/C)^(-siggma)*(Z(+1)/Z)/Pi(+1);
    [name='Definition nominal interest rate), p. 22 top']
    R=1/Q;
    [name='Aggregate output, above eq. (14)']
    Y=A*(N/S)^(1-alppha);
    [name='Definition Real interest rate']
    R=realinterest*Pi(+1);
    @#if money_growth_rule==0
        [name='Monetary Policy Rule, p. 26 bottom/eq. (22)']
        R=1/betta*Pi^phi_pi*(Y/steady_state(Y))^phi_y*exp(nu);
    @#endif
    [name='Market Clearing, eq. (15)']
    C=Y;
    [name='Technology Shock, eq. (6)']
    log(A)=rho_a*log(A(-1))+eps_a;
    [name='Preference Shock, p.54']
    log(Z)=rho_z*log(Z(-1))-eps_z;
    @#if money_growth_rule==0
        [name='Monetary policy shock']
        nu=rho_nu*nu(-1)+eps_nu;
    @#endif
    @#if money_growth_rule==1
        [name='definition nominal money growth']
        money_growth=log(M_real/M_real(-1)*Pi);
        [name='exogenous process for money supply growth rate, eq. (35)']
        money_growth=rho_m*money_growth(-1)+eps_m;
        [name='definition annualized nominal money growth'] 
        money_growth_ann=4*money_growth;
    @#endif
    [name='Definition marginal cost']
    MC=W_real/((1-alppha)*Y/N*S);
    [name='LOM prices, eq. (7)']
    1=theta*Pi^(epsilon-1)+(1-theta)*(Pi_star)^(1-epsilon);
    [name='LOM price dispersion']
    S=(1-theta)*Pi_star^(-epsilon/(1-alppha))+theta*Pi^(epsilon/(1-alppha))*S(-1);
    [name='FOC price setting']
    Pi_star^(1+epsilon*(alppha/(1-alppha)))=x_aux_1/x_aux_2*(1-tau)*epsilon/(epsilon-1);
    [name='Auxiliary price setting recursion 1']
    x_aux_1=Z*C^(-siggma)*Y*MC+betta*theta*Pi(+1)^(epsilon+alppha*epsilon/(1-alppha))*x_aux_1(+1);
    [name='Auxiliary price setting recursion 2']
    x_aux_2=Z*C^(-siggma)*Y+betta*theta*Pi(+1)^(epsilon-1)*x_aux_2(+1);
    [name='Definition log output']
    log_y = log(Y);
    [name='Definition log real wage']
    log_W_real=log(W_real);
    [name='Definition log hours']
    log_N=log(N);
    [name='Annualized inflation']
    pi_ann=4*log(Pi);
    [name='Annualized nominal interest rate']
    i_ann=4*log(R);
    [name='Annualized real interest rate']
    r_real_ann=4*log(realinterest);
    [name='Real money demand, eq. (4)']
    M_real=Y/R^eta;
    [name='definition nominal money stock']
    log_m_nominal=log(M_real*P);
    [name='Definition price level']
    Pi=P/P(-1);
    [name='Definition log price level']
    log_P=log(P);
    [name='Definition log TFP']
    log_A=log(A);
    [name='Definition log preference']
    log_Z=log(Z);
end;

%----------------------------------------------------------------
%  Steady state values
%---------------------------------------------------------------

steady_state_model;
A=1;
Z=1;
S=1;
Pi_star=1;
P=1;
MC=(epsilon-1)/epsilon/(1-tau);
R=1/betta;
Pi=1;
Q=1/R;
realinterest=R;
N=((1-alppha)*MC)^(1/((1-siggma)*alppha+varphi+siggma));
C=A*N^(1-alppha);
W_real=C^siggma*N^varphi;
Y=C;
money_growth=0;
money_growth_ann=0;
nu=0;
x_aux_1=C^(-siggma)*Y*MC/(1-betta*theta*Pi^(epsilon/(1-alppha)));
x_aux_2=C^(-siggma)*Y/(1-betta*theta*Pi^(epsilon-1));
log_y = log(Y);
log_W_real=log(W_real);
log_N=log(N);
pi_ann=4*log(Pi);
i_ann=4*log(R);
r_real_ann=4*log(realinterest);
M_real=Y/R^eta;
log_m_nominal=log(M_real*P);
log_P=log(P);
log_A=0;
log_Z=0;
end;

write_latex_dynamic_model;

resid;
steady;
check;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------

shocks;
    @#if money_growth_rule==0
        var eps_nu = 0.0025^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    @#else   
        var eps_m = 0.0025^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    @#endif
end;

%----------------------------------------------------------------
% generate IRFs for monetary policy shock, replicates Figures 3.1, p. 69 (interest rate rule)
% 3.4, p. 76 (money supply rule)
%----------------------------------------------------------------
@#if money_growth_rule==0
    stoch_simul(order = 1,irf=15) pi_ann log_y log_N log_W_real log_P i_ann r_real_ann log_m_nominal nu;
@#else
    stoch_simul(order = 1,irf=15) pi_ann log_y log_N log_W_real log_P i_ann r_real_ann log_m_nominal money_growth_ann;
@#endif


%----------------------------------------------------------------
% generate IRFs for discount rate shock, replicates Figures 3.2, p. 70 (interest rate rule)
% 3.5, p. 78 (money supply rule)
%----------------------------------------------------------------
shocks;
    @#if money_growth_rule==0
        var eps_nu = 0;   //shut off monetary policy shock
    @#else   
        var eps_m = 0;   //shut off monetary policy shock
    @#endif
var eps_z  = 0.005^2; //unit shock to technology
end;

stoch_simul(order = 1,irf=15,irf_plot_threshold=0) pi_ann log_y log_N log_W_real log_P i_ann r_real_ann log_m_nominal log_Z;


%----------------------------------------------------------------
% generate IRFs, replicates Figures 3.3, p. 73 (interest rate rule)
% 3.6, p. 81 (money supply rule)
%----------------------------------------------------------------
shocks;
var eps_z = 0;   //shut off discount rate shock
var eps_a  = 0.01^2; //unit shock to technology
end;

stoch_simul(order = 1,irf=15,irf_plot_threshold=0) pi_ann log_y log_N log_W_real log_P i_ann r_real_ann log_m_nominal log_A;

% if exist('Gali_2015_chapter_3_results.mat','file')
%     oo_linear=load('Gali_2015_chapter_3_results.mat','oo_')
%     if max(abs(oo_linear.oo_.irfs.y_eps_m/100-oo_.irfs.log_y_eps_m))>1e-10 || ...
%             max(abs(oo_linear.oo_.irfs.y_eps_z/100-oo_.irfs.log_y_eps_z))>1e-10 || ...
%             max(abs(oo_linear.oo_.irfs.y_eps_a/100-oo_.irfs.log_y_eps_a))>1e-10 || ...
%             max(abs(oo_linear.oo_.irfs.p_eps_m/100-oo_.irfs.log_P_eps_m))>1e-10 || ...
%             max(abs(oo_linear.oo_.irfs.p_eps_z/100-oo_.irfs.log_P_eps_z))>1e-10 || ...
%             max(abs(oo_linear.oo_.irfs.p_eps_a/100-oo_.irfs.log_P_eps_a))>1e-10
%         error('IRFs do not match linear model')
%     end
% end
