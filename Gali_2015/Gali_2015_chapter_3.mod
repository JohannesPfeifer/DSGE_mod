/*
 * This file implements the baseline New Keynesian model of Jordi Galí (2015): Monetary Policy, Inflation,
 * and the Business Cycle, Princeton University Press, Second Edition, Chapter 3
 * 
 * THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER
 *
 * Notes:
 *  - all model variables are expressed in deviations from steady state, i.e. in contrast to
 *      to the chapter, both the nominal interest rate and natural output are not in log-levels, but rather mean 0
 *  - in the LOM for the discount rate shock z the shock enters with a minus sign in this mod-file to generate the 
 *      IRF to a -0.5% shock
 *  - the IRF for the nominal rate in Figure 3.6 "Dynamic Responses to a Technology Shock: Money Supply Rule", p. 81
 *      is wrong. It should be identically 0 as can be seen in this mod-file and Galí's slide set accompanying this chapter
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
@#define money_growth_rule=0

var pi          ${\pi}$                 (long_name='inflation')
    y_gap       ${\tilde y}$            (long_name='output gap')
    y_nat       ${y^{nat}}$             (long_name='natural output')      //(in contrast to the textbook defined in deviation from steady state)
    y           ${y}$                   (long_name='output')
    yhat        ${\hat y}$              (long_name='output deviation from steady state')
    r_nat       ${r^{nat}}$             (long_name='natural interest rate')
    r_real      ${r^r}$                 (long_name='real interest rate')     
    i           ${i}$                   (long_name='nominal interest rate')
    n           ${n}$                   (long_name='hours worked')
    m_real      ${m-p}$                 (long_name='real money stock')
    m_growth_ann ${\Delta m}$           (long_name='money growth annualized')
    m_nominal   ${m}$                   (long_name='nominal money stock')
    @#if money_growth_rule==0
        nu      ${\nu}$                 (long_name='AR(1) monetary policy shock process')    
    @#else
        money_growth  ${\Delta m_q}$    (long_name='money growth')
        money_growth_ann  ${\Delta m^{ann}}$    (long_name='money growth annualized')
   @#endif
    a           ${a}$                   (long_name='AR(1) technology shock process')
    r_real_ann  ${r^{r,ann}}$           (long_name='annualized real interest rate')
    i_ann       ${i^{ann}}$             (long_name='annualized nominal interest rate')
    r_nat_ann   ${r^{nat,ann}}$         (long_name='annualized natural interest rate')
    pi_ann      ${\pi^{ann}}$           (long_name='annualized inflation rate')
    z           ${z}$                   (long_name='AR(1) preference shock process')
    p           ${p}$                   (long_name='price level')
    w           ${w}$                   (long_name='nominal wage')
    c           ${c}$                   (long_name='consumption')
    w_real      ${\frac{w}{p}}$         (long_name='real wage')
    mu          ${\mu}$                 (long_name='markup')
    mu_hat      ${\hat \mu}$            (long_name='markup gap')
;     

varexo  eps_a       ${\varepsilon_a}$       (long_name='technology shock')
        @#if money_growth_rule==0
            eps_nu  ${\varepsilon_\nu}$     (long_name='monetary policy shock')
        @#else   
            eps_m       ${\varepsilon_m}$   (long_name='money supply shock innovation')
        @#endif
        eps_z       ${\varepsilon_z}$   (long_name='preference shock innovation')
       ;

parameters alppha       ${\alpha}$     (long_name='capital share')
    betta               ${\beta}$       (long_name='discount factor')
    rho_a               ${\rho_a}$      (long_name='autocorrelation technology shock')
    @#if money_growth_rule==0
        rho_nu          ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')
    @#else   
        rho_m        ${\rho_{\zeta}}$ (long_name='autocorrelation monetary demand shock')
    @#endif
    rho_z               ${\rho_{z}}$    (long_name='autocorrelation preference shock')
    siggma              ${\sigma}$      (long_name='inverse EIS')
    varphi              ${\varphi}$     (long_name='inverse Frisch elasticity')
    phi_pi              ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    phi_y               ${\phi_{y}}$    (long_name='output feedback Taylor Rule')
    eta                 ${\eta}$        (long_name='semi-elasticity of money demand')
    epsilon             ${\epsilon}$    (long_name='demand elasticity')
    theta               ${\theta}$      (long_name='Calvo parameter')
    ;
%----------------------------------------------------------------
% Parametrization, p. 67-75
%----------------------------------------------------------------
siggma = 1;         %p. 67
varphi=5;           %p. 67, Frisch elasticity of 0.2
phi_pi = 1.5;       %p. 68 
phi_y  = 0.125;     %p. 68 (5/4)
theta=3/4;          %p. 67
@#if money_growth_rule==0
    rho_nu =0.5;    %p. 68
@#else   
    rho_m=0.5;      %p. 75
@#endif
rho_z  = 0.5;       %p. 70
rho_a  = 0.9;       %p. 72
betta  = 0.99;      %p. 67
eta  =3.77; %footnote 11, p. 115
alppha=1/4;     	%p. 67
epsilon=9;          %p. 67

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 
//Composite parameters
#Omega=(1-alppha)/(1-alppha+alppha*epsilon);        %defined on page 60
#psi_n_ya=(1+varphi)/(siggma*(1-alppha)+varphi+alppha);   %defined on page 62
#lambda=(1-theta)*(1-betta*theta)/theta*Omega;      %defined on page 61
#kappa=lambda*(siggma+(varphi+alppha)/(1-alppha));     %defined on page 63
[name='New Keynesian Phillips Curve eq. (22)']
pi=betta*pi(+1)+kappa*y_gap;
[name='Dynamic IS Curve eq. (23)']
y_gap=-1/siggma*(i-pi(+1)-r_nat)+y_gap(+1);
@#if money_growth_rule==0
    [name='Interest Rate Rule eq. (26)']
    i=phi_pi*pi+phi_y*yhat+nu;
@#endif
[name='Definition natural rate of interest eq. (24)']
r_nat=-siggma*psi_n_ya*(1-rho_a)*a+(1-rho_z)*z;
[name='Definition real interest rate']
r_real=i-pi(+1);
[name='Definition natural output, eq. (20)']
y_nat=psi_n_ya*a;
[name='Definition output gap']
y_gap=y-y_nat;
@#if money_growth_rule==0
    [name='Monetary policy shock']
    nu=rho_nu*nu(-1)+eps_nu;
@#endif
[name='TFP shock']
a=rho_a*a(-1)+eps_a;
[name='Production function (eq. 14)']
y=a+(1-alppha)*n;
[name='Preference shock, p. 54']
z     = rho_z*z(-1) - eps_z;
[name='Money growth (derived from eq. (4))']
m_growth_ann=4*(y-y(-1)-eta*(i-i(-1))+pi);
@#if money_growth_rule==1
    [name='Real money demand, eq. (4)']
    m_real=y-eta*i;
    [name='definition nominal money growth']
    money_growth=m_real-m_real(-1)+pi;
    [name='exogenous process for money supply growth rate, eq. (35)']
    money_growth=rho_m*money_growth(-1)+eps_m;
    [name='definition annualized nominal money growth'] 
    money_growth_ann=4*money_growth;
@#else
    [name='Real money demand (eq. 4)']
    m_real=y-eta*i;
@#endif
[name='Annualized nominal interest rate']
i_ann=4*i;
[name='Annualized real interest rate']
r_real_ann=4*r_real;
[name='Annualized natural interest rate']
r_nat_ann=4*r_nat;
[name='Annualized inflation']
pi_ann=4*pi;
[name='Output deviation from steady state']
yhat=y-steady_state(y);
[name='Definition price level']
pi=p-p(-1);
[name='resource constraint, eq. (12)']
y=c;
[name='FOC labor, eq. (2)']
w-p=siggma*c+varphi*n;
[name='definition real wage']
w_real=w-p;
[name='definition nominal money stock']
m_nominal=m_real+p;
[name='average price markup, eq. (18)']
mu=-(siggma+(varphi+alppha)/(1-alppha))*y+(1+varphi)/(1-alppha)*a;
[name='average price markup, eq. (20)']
mu_hat=-(siggma+(varphi+alppha)/(1-alppha))*y_gap;
end;


%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------


shocks;
    @#if money_growth_rule==0
        var eps_nu = 0.25^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    @#else   
        var eps_m = 0.25^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    @#endif
end;

%----------------------------------------------------------------
%  steady states: all 0 due to linear model
%---------------------------------------------------------------
resid;
steady;
check;

%----------------------------------------------------------------
% generate IRFs for monetary policy shock, replicates Figures 3.1, p. 69 (interest rate rule)
% 3.4, p. 76 (money supply rule)
%----------------------------------------------------------------
@#if money_growth_rule==0
    stoch_simul(order = 1,irf=15) y_gap pi_ann y n w_real p i_ann r_real_ann m_nominal nu;
@#else
    stoch_simul(order = 1,irf=15) y_gap pi_ann y n w_real p i_ann r_real_ann m_nominal money_growth_ann;
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
var eps_z  = 0.5^2; //unit shock to preferences 
end;

stoch_simul(order = 1,irf=15,irf_plot_threshold=0) y_gap pi_ann y n w_real p i_ann r_real_ann m_nominal z;


%----------------------------------------------------------------
% generate IRFs, replicates Figures 3.3, p. 73 (interest rate rule)
% 3.6, p. 81 (money supply rule)
%----------------------------------------------------------------
shocks;
    @#if money_growth_rule==0
        var eps_z = 0;   //shut off preference shock
    @#else   
        var eps_z = 0;   //shut off preference shock
    @#endif
var eps_a  = 1^2; //unit shock to technology
end;

stoch_simul(order = 1,irf=15,irf_plot_threshold=0) y_gap pi_ann y n w_real p i_ann r_real_ann m_nominal a;

