/*
 * This file implements the optimal monetary policy at the ZLB under discretion exercise 
 * in Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, 
 * Princeton University Press, Second Edition, Chapter 5.4.1
 *
 * Notes:
 *  - all model variables are expressed in deviations from steady state, except for 
 *      the interest rates.
 *  - This file replicates the IRFs for discretion reported in Figure 5.3
 *  - The dummy-variable ZLB_indicator is used to shut off the Taylor rule during the 
 *      ZLB period and to switch it on again after the ZLB-period. This is necessary,
 *      because without a Taylor rule to govern monetary policy after the ZLB period,
 *      there is scope for multiple equilibria. Switching on the Taylor rule assures that
 *      the solution pi=x=0 is selected. Otherwise, one would need to set the number of 
 *      simulation periods to 6 so that the terminal condition of 
 *      all variables being back to steady state is imposed immediately after 
 *      the end of the ZLB period. 
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
    x           ${x}$           (long_name='welfare-relevant output gap')
    i           ${i}$           (long_name='nominal interest rate')
    r_nat_ann   ${r^{nat,ann}}$ (long_name='annualized natural interest rate')
    pi_ann      ${\pi^{ann}}$   (long_name='annualized inflation rate')
    p           ${p}$           (long_name='price level')
    xi_2        ${\xi_2}$       (long_name='Langrange multiplier')
    i_ann       ${i^{ann}}$     (long_name='annualized nominal interest rate')  
    i_taylor    ${i^{Taylor}}$  (long_name='nominal interest rate Taylor rule')  
    ;     

varexo r_nat ${r^n}$   (long_name='natural rate of interest')
        ZLB_indicator   ${ZLB\_IND}$ (long_name='ZLB Indicator');

parameters alppha   ${\alppha}$         (long_name='capital share')
        betta       ${\beta}$           (long_name='discount factor')
        siggma      ${\sigma}$          (long_name='log utility')
        varphi      ${\varphi}$         (long_name='unitary Frisch elasticity')
        epsilon     ${\epsilon}$        (long_name='demand elasticity')
        theta       ${\theta}$          (long_name='Calvo parameter')
        phi_pi      ${\phi_\pi}$        (long_name='Taylor rule inflation feedback')
    ;
%----------------------------------------------------------------
% Parametrization, p. 52
%----------------------------------------------------------------
siggma = 1;
varphi=5;
theta=3/4;
betta  = 0.99;
alppha = 1/4;
epsilon= 9;
phi_pi=1.5;

model; 
#Omega=(1-alppha)/(1-alppha+alppha*epsilon);         %defined on page 60
#lambda=(1-theta)*(1-betta*theta)/theta*Omega;       %defined on page 61
#kappa=lambda*(siggma+(varphi+alppha)/(1-alppha));   %defined on page 63
#vartheta=kappa/epsilon;                             %defined on page 128
[name='New Keynesian Phillips Curve eq. (29)']
pi=betta*pi(+1)+kappa*x;
[name='Dynamic IS Curve eq. (30)']
x=x(+1)-1/siggma*(max(i*ZLB_indicator+(1-ZLB_indicator)*i_taylor,0)-pi(+1)-r_nat);
[name='FOC, eq. (33)']
vartheta*x=-kappa*pi-xi_2;
[name='Complementary slackness condition']
max(xi_2,0)*max(i*ZLB_indicator+(1-ZLB_indicator)*i_taylor,0)=0;
[name='Taylor rule']
i_taylor=1+phi_pi*pi;
[name='Annualized natural interest rate']
r_nat_ann=4*r_nat;
[name='Annualized inflation']
pi_ann=4*pi;
[name='Definition price level']
pi=p-p(-1);
[name='Annualized nominal interest rate']
i_ann=4*max(i*ZLB_indicator+(1-ZLB_indicator)*i_taylor,0);
end;

%%set initial and terminal condition to steady state of non-ZLB period
initval;
r_nat=1;
i=1;
i_taylor=1;
i_ann=4*i;
r_nat_ann=4*r_nat;
ZLB_indicator=0;
end;

steady;

%define natural rate shock
shocks;
var r_nat;
periods 1:6;
values -1;
var ZLB_indicator;
periods 1:6;
values 1;
end;

perfect_foresight_setup(periods=20);
perfect_foresight_solver(stack_solve_algo=7,solve_algo=1);

options_.rplottype=2;
rplot x pi_ann i_ann r_nat_ann;