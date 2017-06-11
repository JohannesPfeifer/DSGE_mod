/*
 * This file implements the New Keynesian model with price and wage rigidities under optimal policy 
 * with commitment (Ramsey) of Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton 
 * University Press, Second Edition, Chapter 6.4
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER
 *
 * Notes:
 *  - all model variables are expressed in deviations from steady state, i.e. in contrast to
 *      to the chapter, the nominal interest rate, natural output, and the natural real wage are not in log-levels, but rather mean 0
 *  - in the LOM for the discount rate shock z the shock enters with a minus sign in this mod-file to generate the 
 *      IRF to a -0.5% shock
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

var pi_p        ${\pi^p}$               (long_name='price inflation')
    y_gap       ${\tilde y}$            (long_name='output gap')
    y_nat       ${y^{nat}}$             (long_name='natural output')      //(in contrast to the textbook defined in deviation from steady state)
    y           ${y}$                   (long_name='output')
    yhat        ${\hat y}$              (long_name='output deviation from steady state')
    r_nat       ${r^{nat}}$             (long_name='natural interest rate')
    r_real      ${r^r}$                 (long_name='real interest rate')     
    i           ${i}$                   (long_name='nominal interrst rate')
    n           ${n}$                   (long_name='hours worked')
    m_real      ${(m-p)}$                 (long_name='real money stock')
    m_growth_ann ${\Delta m}$           (long_name='money growth annualized')
    m_nominal   ${m}$                   (long_name='nominal money stock')
    nu          ${\nu}$                 (long_name='AR(1) monetary policy shock process')    
    a           ${a}$                   (long_name='AR(1) technology shock process')
    r_real_ann  ${r^{r,ann}}$           (long_name='annualized real interest rate')
    i_ann       ${i^{ann}}$             (long_name='annualized nominal interest rate')
    r_nat_ann   ${r^{nat,ann}}$         (long_name='annualized natural interest rate')
    pi_p_ann    ${\pi^{p,ann}}$         (long_name='annualized inflation rate')
    z           ${z}$                   (long_name='AR(1) preference shock process')
    p           ${p}$                   (long_name='price level')
    w           ${w}$                   (long_name='nominal wage')
    c           ${c}$                   (long_name='consumption')
    w_real      $\omega$                (long_name='real wage')
    w_gap       ${\tilde \omega}$       (long_name='real wage gap')
    pi_w        ${\pi^w}$               (long_name='wage inflation')
    w_nat       ${w^{nat}}$             (long_name='natural real wage')
    mu_p        ${\mu^p}$               (long_name='markup')
    pi_w_ann    ${\pi^{w,ann}}$         (long_name='annualized wage inflation rate')
;     

varexo  eps_a       ${\varepsilon_a}$   (long_name='technology shock')
        eps_nu      ${\varepsilon_\nu}$ (long_name='monetary policy shock')
        eps_z       ${\varepsilon_z}$   (long_name='preference shock innovation')
       ;

parameters alppha       ${\alpha}$     (long_name='capital share')
    betta               ${\beta}$       (long_name='discount factor')
    rho_a               ${\rho_a}$      (long_name='autocorrelation technology shock')
    rho_nu              ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')
    rho_z               ${\rho_{z}}$    (long_name='autocorrelation monetary demand shock')
    siggma              ${\sigma}$      (long_name='inverse EIS')
    varphi              ${\varphi}$     (long_name='inverse Frisch elasticity')
    phi_pi              ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    phi_y               ${\phi_{y}}$    (long_name='output feedback Taylor Rule')
    eta                 ${\eta}$        (long_name='semi-elasticity of money demand')
    epsilon_p           ${\epsilon_p}$  (long_name='demand elasticity goods')
    theta_p             ${\theta_p}$    (long_name='Calvo parameter prices')
    epsilon_w           ${\epsilon_w}$  (long_name='demand elasticity labor services')
    theta_w             ${\theta_w}$    (long_name='Calvo parameter wages')
    lambda_p            ${\lambda_p}$   (long_name='composite parameter Phillips Curve')
    lambda_w            ${\lambda_w}$   (long_name='composite parameter wage Phillips Curve')
    ;
%----------------------------------------------------------------
% Parametrization, p. 67  and p. 113-115
%----------------------------------------------------------------
siggma = 1;
varphi=5;
phi_pi = 1.5;
phi_y  = 0.125;
theta_p=3/4;
rho_nu =0.5;
rho_z  = 0.5;
rho_a  = 0.9;
betta  = 0.99;
eta  =3.77; %footnote 11, p. 115
alppha=1/4;
epsilon_p=9;

epsilon_w=4.5;
theta_w=3/4;
%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 
//Composite parameters
#Omega=(1-alppha)/(1-alppha+alppha*epsilon_p);              %defined on page 166
#psi_n_ya=(1+varphi)/(siggma*(1-alppha)+varphi+alppha);     %defined on page 171
#psi_n_wa=(1-alppha*psi_n_ya)/(1-alppha);                   %defined on page 171
#aleph_p=alppha*lambda_p/(1-alppha);                         %defined on page 172
#aleph_w=lambda_w*(siggma+varphi/(1-alppha));               %defined on page 172
[name='New Keynesian Phillips Curve eq. (18)']
pi_p=betta*pi_p(+1)+aleph_p*y_gap+lambda_p*w_gap;
[name='New Keynesian Wage Phillips Curve eq. (22)']
pi_w=betta*pi_w(+1)+aleph_w*y_gap-lambda_w*w_gap;
[name='Dynamic IS Curve eq. (22)']
y_gap=-1/siggma*(i-pi_p(+1)-r_nat)+y_gap(+1);
[name='Definition natural rate of interest eq. (24)']
r_nat=-siggma*psi_n_ya*(1-rho_a)*a+(1-rho_z)*z;
w_gap=w_gap(-1)+pi_w-pi_p-(w_nat-w_nat(-1));
[name='Definition natural wage, eq (16)']
w_nat=psi_n_wa*a;
[name='Definition markup']
mu_p=-alppha/(1-alppha)*y_gap-w_gap;
[name='Definition real wage gap, p. 171']
w_gap=w_real-w_nat;
[name='Definition real interest rate']
r_real=i-pi_p(+1);
[name='Definition natural output, eq. (20)']
y_nat=psi_n_ya*a;
[name='Definition output gap']
y_gap=y-y_nat;
[name='Monetary policy shock']
nu=rho_nu*nu(-1)+eps_nu;
[name='TFP shock']
a=rho_a*a(-1)+eps_a;
[name='Production function, p. 171']
y=a+(1-alppha)*n;
[name='Preference shock, p. 54']
z     = rho_z*z(-1) - eps_z;
[name='Money growth (derived from eq. (4))']
m_growth_ann=4*(y-y(-1)-eta*(i-i(-1))+pi_p);
[name='Real money demand (eq. 4)']
m_real=y-eta*i;
[name='Annualized nominal interest rate']
i_ann=4*i;
[name='Annualized real interest rate']
r_real_ann=4*r_real;
[name='Annualized natural interest rate']
r_nat_ann=4*r_nat;
[name='Annualized inflation']
pi_p_ann=4*pi_p;
[name='Annualized wage inflation']
pi_w_ann=4*pi_w;
[name='Output deviation from steady state']
yhat=y-steady_state(y);
[name='Definition price level']
pi_p=p-p(-1);
[name='resource constraint, eq. (12)']
y=c;
[name='definition real wage']
w_real=w-p;
[name='definition real wage']
m_nominal=m_real+p;


end;

steady_state_model;
lambda_p=(1-theta_p)*(1-betta*theta_p)/theta_p*((1-alppha)/(1-alppha+alppha*epsilon_p));      %defined on page 166
lambda_w=(1-theta_w)*(1-betta*theta_w)/(theta_w*(1+epsilon_w*varphi));      %defined on page 170
end;
%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------


shocks;
    var eps_a       = 1; 
    var eps_z       = 1; 
end;


%----------------------------------------------------------------
% generate IRFs for technology shock under optimal policy, replicates Figures 6.3, p. 182
%----------------------------------------------------------------
//planner objective, uses lambda_w and lambda_p updated in steady_state_model-block
planner_objective 0.5*((siggma+(varphi+alppha)/(1-alppha))*y_gap^2+ epsilon_p/lambda_p*pi_p^2+epsilon_w*(1-alppha)/lambda_w*pi_w^2);

ramsey_policy(instruments=(i),irf=16,planner_discount=betta,noprint) y_gap pi_p_ann pi_w_ann w_real;
oo_baseline=oo_;

%flexible wage case
set_param_value('theta_w',0.0000000001);
set_param_value('theta_p',3/4);
ramsey_policy(instruments=(i),irf=16,planner_discount=betta,noprint) y_gap pi_p_ann pi_w_ann w_real;
oo_flexible_wages=oo_;

%flexible price case
set_param_value('theta_w',3/4)
set_param_value('theta_p',0.000000001)
ramsey_policy(instruments=(i),irf=16,planner_discount=betta,noprint) y_gap pi_p_ann pi_w_ann w_real;
oo_flexible_prices=oo_;


figure('Name','Dynamic Responses to a technology shock under optimal policy')
subplot(2,2,1)
plot(1:options_.irf,oo_baseline.irfs.y_gap_eps_a,'-o',1:options_.irf,oo_flexible_wages.irfs.y_gap_eps_a,'-d',1:options_.irf,oo_flexible_prices.irfs.y_gap_eps_a,'-s')
ylim([-0.1 0.1])
title('Output gap')
ll=legend('baseline','flexible wages','flexible prices');
set(ll,'Location','SouthEast');
subplot(2,2,2)
plot(1:options_.irf,oo_baseline.irfs.pi_p_ann_eps_a,'-o',1:options_.irf,oo_flexible_wages.irfs.pi_p_ann_eps_a,'-d',1:options_.irf,oo_flexible_prices.irfs.pi_p_ann_eps_a,'-s')
title('Price Inflation')
subplot(2,2,3)
plot(1:options_.irf,oo_baseline.irfs.pi_w_ann_eps_a,'-o',1:options_.irf,oo_flexible_wages.irfs.pi_w_ann_eps_a,'-d',1:options_.irf,oo_flexible_prices.irfs.pi_w_ann_eps_a,'-s')
title('Wage inflation')
subplot(2,2,4)
plot(1:options_.irf,oo_baseline.irfs.w_real_eps_a,'-o',1:options_.irf,oo_flexible_wages.irfs.w_real_eps_a,'-d',1:options_.irf,oo_flexible_prices.irfs.w_real_eps_a,'-s')
title('Real wage')


%----------------------------------------------------------------
% generate first row of Table 6.1, p. 186
%----------------------------------------------------------------
shocks;
    var eps_a       = 1; 
end;
set_param_value('theta_w',3/4);
set_param_value('theta_p',3/4);

ramsey_policy(instruments=(i),irf=16,planner_discount=betta) y_gap pi_p pi_w;
oo_baseline=oo_;

y_gap_pos=strmatch('y_gap',var_list_ ,'exact');
pi_p_pos=strmatch('pi_p',var_list_ ,'exact');
pi_w_pos=strmatch('pi_w',var_list_ ,'exact');

%read out current parameter values
par.alppha=M_.params(strmatch('alppha',M_.param_names,'exact'));
par.epsilon_p=M_.params(strmatch('epsilon_p',M_.param_names,'exact'));
par.epsilon_w=M_.params(strmatch('epsilon_w',M_.param_names,'exact'));
par.siggma=M_.params(strmatch('siggma',M_.param_names,'exact'));
par.varphi=M_.params(strmatch('varphi',M_.param_names,'exact'));
par.lambda_w=M_.params(strmatch('lambda_w',M_.param_names,'exact'));
par.lambda_p=M_.params(strmatch('lambda_p',M_.param_names,'exact'));

variance.y_gap=oo_.var(y_gap_pos,y_gap_pos);
variance.pi_p=oo_.var(pi_p_pos,pi_p_pos);
variance.pi_w=oo_.var(pi_w_pos,pi_w_pos);
L=0.5*((par.siggma+(par.varphi+par.alppha)/(1-par.alppha))*variance.y_gap+ par.epsilon_p/par.lambda_p*variance.pi_p+par.epsilon_w*(1-par.alppha)/par.lambda_w*variance.pi_w)

labels=strvcat('sigma(pi_p)','sigma(pi_w)','sigma(tilde y)','L');
headers=strvcat(' ','Optimal');
values=[sqrt([variance.pi_p;variance.pi_w;variance.y_gap]);L];
options_.noprint=0;
dyntable(options_,'Evaluation of Simple Rules',headers,labels,values,size(labels,2)+2,4,3)
