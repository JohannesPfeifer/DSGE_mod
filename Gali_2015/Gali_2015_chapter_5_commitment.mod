/*
 * This file implements the optimal monetary policy under commitment exercise 
 * in Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, 
 * Princeton University Press, Second Edition, Chapter 5.2.2
 *
 * It demonstrates how to use the ramsey_policy command of Dynare.
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER
 *
 * Notes:
 *      - all model variables are expressed in deviations from steady state, i.e. 
 *        in contrast to to the chapter, both the nominal interest rate and 
 *        natural output are not in log-levels, but rather mean 0
 *      - ramsey_policy computes the optimal interest rate under commitment
 *      - This file replicates the IRFs reported in Figures 5.1 and 5.2
 *      - An undistorted steady state is assumed in this chapter
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


var pi          ${\pi}$         (long_name='inflation')
    y_gap       ${tilde y}$     (long_name='output gap')
    y_nat       ${y^{nat}}$     (long_name='natural output')
    y           ${y}$           (long_name='output')
    r_e         ${r^{e}}$       (long_name='efficient interest rate')
    y_e         ${y^{nat}}$     (long_name='efficient output') 
    x           ${x}$           (long_name='welfare-relevant output gap')
    r_real      ${r^r}$         (long_name='real interest rate')     
    i           ${i}$           (long_name='nominal interest rate')
    n           ${n}$           (long_name='hours worked')
    m_growth_ann ${\Delta m}$   (long_name='money growth')
    u           ${u}$           (long_name='AR(1) cost push shock process')
    a           ${a}$           (long_name='AR(1) technology shock process')
    r_real_ann  ${r^{r,ann}}$   (long_name='annualized real interest rate')
    i_ann       ${i^{ann}}$     (long_name='annualized nominal interest rate')
    pi_ann      ${\pi^{ann}}$   (long_name='annualized inflation rate')
    p           ${p}$           (long_name='price level')
    z           ${z}$           (long_name='AR(1) preference shock process')
    ;     

varexo eps_a ${\varepsilon_a}$   (long_name='technology shock')
       eps_u ${\varepsilon_u}$   (long_name='monetary policy shock')
       eps_z ${\varepsilon_z}$   (long_name='preference shock innovation')
;

parameters alppha   ${\alppha}$         (long_name='capital share')
        betta       ${\beta}$           (long_name='discount factor')
        rho_a       ${\rho_a}$          (long_name='autocorrelation technology shock')
        rho_u       ${\rho_{u}}$        (long_name='autocorrelation cost push shock')
        rho_z       ${\rho_{z}}$        (long_name='autocorrelation monetary demand shock')
        siggma      ${\sigma}$          (long_name='log utility')
        varphi      ${\varphi}$         (long_name='unitary Frisch elasticity')
        eta         ${\eta}$            (long_name='semi-elasticity of money demand')
        epsilon     ${\epsilon}$        (long_name='demand elasticity')
        theta       ${\theta}$          (long_name='Calvo parameter')
        Omega       ${\Omega}$          (long_name='Composite parameter Phillips curve')
        lambda      ${\lambda}$         (long_name='Composite parameter Phillips curve')
        kappa       ${\kappa}$          (long_name='Composite parameter Phillips curve')
        vartheta    ${\vartheta}$       (long_name='weight of x in utility function')
    ;
%----------------------------------------------------------------
% Parametrization, p. 52
%----------------------------------------------------------------
rho_u=0;
siggma = 1;
varphi=5;
theta=3/4;
rho_z  = 0.5;
rho_a  = 0.9;
betta  = 0.99;
eta    = 3.77; %footnote 11, p. 115
alppha = 1/4;
epsilon= 9;



%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 
//Composite parameters
#psi_n_ya=(1+varphi)/(siggma*(1-alppha)+varphi+alppha); %defined on page 62
[name='Definition efficient interest rate, below equation (7)']
r_e=siggma*(y_e(+1)-y_e)+(1-rho_z)*z;
[name='Definition efficient output']
y_e=psi_n_ya*a;
[name='Definition linking various output gaps, middle page 128']
y_gap=x+(y_e-y_nat);
[name='New Keynesian Phillips Curve eq. (2)']
pi=betta*pi(+1)+kappa*x + u;
[name='Dynamic IS Curve eq. (7)']
x=x(+1)-1/siggma*(i-pi(+1)-r_e);
[name='Definition real interest rate']
r_real=i-pi(+1);
[name='Implicit definition of natural output, following from definition of u']
u=kappa*(y_e-y_nat);
[name='Definition output gap']
y_gap=y-y_nat;
[name='cost push shock, equation (3)']
u=rho_u*u(-1)+eps_u;
[name='TFP shock']
a=rho_a*a(-1)+eps_a;
[name='Production function, Ch. 3, (eq. 14)']
y=a+(1-alppha)*n;
[name='Money growth (derived from Ch. 3, eq. (4))']
m_growth_ann=4*(y-y(-1)-eta*(i-i(-1))+pi);
[name='Annualized nominal interest rate']
i_ann=4*i;
[name='Annualized real interest rate']
r_real_ann=4*r_real;
[name='Annualized inflation']
pi_ann=4*pi;
[name='Definition price level']
pi=p-p(-1);
[name='Preference shock, p. 54']
z     = rho_z*z(-1) - eps_z;

end;

steady_state_model;
Omega=(1-alppha)/(1-alppha+alppha*epsilon);         %defined on page 60
lambda=(1-theta)*(1-betta*theta)/theta*Omega;       %defined on page 61
kappa=lambda*(siggma+(varphi+alppha)/(1-alppha));   %defined on page 63
vartheta=kappa/epsilon;                             %defined on page 128
end;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------

shocks;
var eps_u = 1;
// var eps_z = 1;
end;

%----------------------------------------------------------------
%  Replicate Figure 5.1 under commitment: response to 
%  transitory cost-push shock
%---------------------------------------------------------------

planner_objective pi^2 +vartheta*x^2;
ramsey_policy(instruments=(i),irf=13,planner_discount=betta) x pi p u;

%----------------------------------------------------------------
%  Replicate Figure 5.2 under commitment: response to 
%  persistent cost-push shock
%---------------------------------------------------------------

set_param_value('rho_u',0.8);
ramsey_policy(instruments=(i),irf=13) x pi p u;
