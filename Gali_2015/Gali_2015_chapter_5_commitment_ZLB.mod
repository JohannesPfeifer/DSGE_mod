/*
 * This file implements the optimal monetary policy at the ZLB under commitment exercise 
 * in Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, 
 * Princeton University Press, Second Edition, Chapter 5.4.2
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5 (I.E. THE CURRENT UNSTABLE VERSION)
 *
 * Notes:
 *  - This mod-file makes use of the Levenberg-Marquardt mixed complementarity problem (lmmcp)
 *      approach to deal with the ZLB constraint. For this purpose, the Lagrangian is written to 
 *      include the instrument i (in contrast to Gali's setup where he uses a inequality constraint obtained 
 *      by substituting out i). I.e. we have
 *      \[
 *          L = \sum\limits_{t = 0}^\infty \beta ^t}\left[ \frac{1}{2}\left(\pi_t^2 + \theta x_t^2 \right) + \xi_{1,t}\left(\pi_t - \kappa x_t - \beta \pi_{t+1} \right) 
 *              + \xi_{2,t}\left(x_t - x_{t+1} + \frac{1}{\sigma}\left(i_t - \pi_{t+1} - r_t^n} \right) \right) \right] 
 *      \]
 *      which gives rise to the FOC
 *      \[
 *          \xi_{2,t}\frac{1}{\sigma} = 0
 *      \]
 *      This FOC has to hold whenever the non-negativity constraint i>0 is not binding. This fact is communicated to
 *      the LMMCP solver by attaching the mcp-tag mcp='i>0' to the equation. Thereby we avoid using the complementary
 *      slackness condition $\xi_{2,t}*i=0$ that would give rise to a singular Jacobian.
 *  - all model variables are expressed in deviations from steady state, except for 
 *      the interest rates.
 *  - This file replicates the IRFs for commitment reported in Figure 5.3
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
    xi_1        ${\xi_1}$       (long_name='Langrange multiplier')
    xi_2        ${\xi_2}$       (long_name='Langrange multiplier')
    i_ann       ${i^{ann}}$     (long_name='annualized nominal interest rate')
    ;     

varexo r_nat ${r^n}$   (long_name='natural rate of interest');

parameters alppha   ${\alppha}$         (long_name='capital share')
        betta       ${\beta}$           (long_name='discount factor')
        siggma      ${\sigma}$          (long_name='log utility')
        varphi      ${\varphi}$         (long_name='unitary Frisch elasticity')
        epsilon     ${\epsilon}$        (long_name='demand elasticity')
        theta       ${\theta}$          (long_name='Calvo parameter')
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


model; 
#Omega=(1-alppha)/(1-alppha+alppha*epsilon);         %defined on page 60
#lambda=(1-theta)*(1-betta*theta)/theta*Omega;       %defined on page 61
#kappa=lambda*(siggma+(varphi+alppha)/(1-alppha));   %defined on page 63
#vartheta=kappa/epsilon;                             %defined on page 128
[name='New Keynesian Phillips Curve eq. (29)']
pi=betta*pi(+1)+kappa*x;
[name='Dynamic IS Curve eq. (30)']
x=x(+1)-1/siggma*(i-pi(+1)-r_nat);
[name='FOC 2, eq. (35)']
pi + xi_1 - xi_1(-1)-1/(betta*siggma)*xi_2(-1)=0;
[name='FOC 2, eq. (36)']
vartheta*x-kappa*xi_1+xi_2-1/betta*xi_2(-1)=0;
[name='FOC w.r.t. to i',mcp='i>0']
xi_2*1/siggma=0;
[name='Annualized natural interest rate']
r_nat_ann=4*r_nat;
[name='Annualized inflation']
pi_ann=4*pi;
[name='Definition price level']
pi=p-p(-1);
[name='Annualized nominal interest rate']
i_ann=4*max(i,0);
end;

%%set initial and terminal condition to steady state of non-ZLB period
initval;
r_nat=1;
i=1;
i_ann=4*i;
r_nat_ann=4*r_nat;
end;
steady;

%define natural rate shock
shocks;
var r_nat;
periods 1:6;
values -1;
end;

perfect_foresight_setup(periods=50);
perfect_foresight_solver(lmmcp);


figure
subplot(2,2,1)
plot(0:12,oo_.endo_simul(strmatch('x',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+13),'-x')
axis([0 12 -15 5])

subplot(2,2,2)
plot(0:12,oo_.endo_simul(strmatch('pi',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+13),'-x')
axis([0 12 -25 5])

subplot(2,2,3)
plot(0:12,oo_.endo_simul(strmatch('i_ann',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+13),'-x')
axis([0 12 -2 6])

subplot(2,2,4)
plot(0:12,oo_.endo_simul(strmatch('r_nat_ann',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+13),'-')
axis([0 12 -6 6])