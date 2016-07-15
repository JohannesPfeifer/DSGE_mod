/*
 * This file implements the optimal monetary policy under discretion exercise 
 * in Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, 
 * Princeton University Press, Second Edition, Chapter 5.2.1
 *
 * It demonstrates how to use the discretionary_policy command of Dynare.
 *
 *  THIS MOD-FILE REQUIRES DYNARE 4.5 (I.E. THE CURRENT UNSTABLE VERSION)
 *
 * Notes:
 *      - all model variables are expressed in deviations from steady state, i.e. 
 *        in contrast to to the chapter, both the nominal interest rate and 
 *        natural output are not in log-levels, but rather mean 0
 *      - discretionary policy computes the optimal interest rate under discretion; 
 *        the mod-file has the optimal analytical rule commented out as the last 
 *        equation. Commenting it in allows using stoch_simul to generate the same
 *        IRFs
 *      - This file replicates the IRFs reported in Figures 5.1 and 5.2
 *      - Note that the output response under discretion in Figures 5.1 and 5.2 in the 
 *          book is wrong, while the one computed by the mod-file is correct. This can be
 *          easily verified by using equation (5) on page 130:
 *              x=-kappa/(kappa^2+vartheta*(1-betta*rho_u))*u
 *          and is done below
 *      - Section 5.2 has an undistorted steady state so that the efficient and the 
 *        natural level of output coincide.
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
    r_nat       ${r^{nat}}$     (long_name='natural interest rate')
    r_real      ${r^r}$         (long_name='real interest rate')     
    i           ${i}$           (long_name='nominal interest rate')
    n           ${n}$           (long_name='hours worked')
    m_growth_ann ${\Delta m}$   (long_name='money growth')
    u           ${u}$           (long_name='AR(1) cost push shock process')
    a           ${a}$           (long_name='AR(1) technology shock process')
    r_real_ann  ${r^{r,ann}}$   (long_name='annualized real interest rate')
    i_ann       ${i^{ann}}$     (long_name='annualized nominal interest rate')
    r_nat_ann   ${r^{nat,ann}}$ (long_name='annualized natural interest rate')
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
        phi_pi      ${\phi_\pi}$        (long_name='Calvo parameter')
        Omega       ${\Omega}$          (long_name='Composite parameter Phillips curve')
        lambda      ${\lambda}$         (long_name='Composite parameter Phillips curve')
        kappa       ${\kappa}$          (long_name='Composite parameter Phillips curve')
        vartheta    ${\vartheta}$       (long_name='weight of x in utility function')
        Theta_i     ${\Theta_i}$        (long_name='weight of u in optimal monetary policy rule')
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
phi_pi=1.5;


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
[name='Definition natural rate of interest, Ch. 3, eq. (24)']
r_nat=-siggma*psi_n_ya*(1-rho_a)*a+(1-rho_z)*z;
[name='Definition real interest rate']
r_real=i-pi(+1);
[name='Definition natural output, eq. (20)']
y_nat=psi_n_ya*a;
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
[name='Annualized natural interest rate']
r_nat_ann=4*r_nat;
[name='Annualized inflation']
pi_ann=4*pi;
[name='Definition price level']
pi=p-p(-1);
[name='Preference shock, p. 54']
z     = rho_z*z(-1) - eps_z;
% [name='Interest Rate Rule that implements optimal solution, eq. (10)']
% i=r_e+phi_pi*pi+Theta_i*u;
end;

steady_state_model;
Omega=(1-alppha)/(1-alppha+alppha*epsilon);         %defined on page 60
lambda=(1-theta)*(1-betta*theta)/theta*Omega;       %defined on page 61
kappa=lambda*(siggma+(varphi+alppha)/(1-alppha));   %defined on page 63
vartheta=kappa/epsilon;                             %defined on page 128
Theta_i=(siggma*kappa*(1-rho_u)-vartheta*(phi_pi-rho_u))/(kappa^2+vartheta*(1-betta*rho_u)); %defined on page 133
end;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------

shocks;
var eps_u = 1;
// var eps_z = 1;
end;

//planner objective with parameters updated in steady state file
planner_objective pi^2 +vartheta*x^2;
discretionary_policy(instruments=(i),irf=13,planner_discount=betta,discretionary_tol=1e-12) x pi p u;

%Verify result using analytical solution
par.kappa=M_.params(strmatch('kappa',M_.param_names,'exact'));
par.vartheta=M_.params(strmatch('vartheta',M_.param_names,'exact'));
par.betta=M_.params(strmatch('betta',M_.param_names,'exact'));
par.rho_u=M_.params(strmatch('rho_u',M_.param_names,'exact'));
x_predicted=-par.kappa/(par.kappa^2+par.vartheta*(1-par.betta*par.rho_u));
fprintf('Analytical output response: %10.8f\n',x_predicted)
fprintf('Computed output response: %10.8f\n',oo_.irfs.x_eps_u(1))
figure
subplot(2,2,1)
plot(0:options_.irf-1,oo_.irfs.x_eps_u,'-o')
axis tight
title('Output gap')
legend('Discretion','Location','South')
subplot(2,2,2)
plot(0:options_.irf-1,oo_.irfs.pi_eps_u,'-o')
axis tight
title('Inflation')
subplot(2,2,3)
plot(0:options_.irf-1,oo_.irfs.p_eps_u,'-o')
axis tight
ylim([-0.1 0.5])
title('Price level')
subplot(2,2,4)
plot(0:options_.irf-1,oo_.irfs.u_eps_u,'-o')
axis tight
title('Cost-push shock')
print -depsc2 Figure_5_1_discretion

set_param_value('rho_u',0.8);
discretionary_policy(instruments=(i),irf=13,discretionary_tol=1e-12) x pi p u;
%Verify result using analytical solution
par.kappa=M_.params(strmatch('kappa',M_.param_names,'exact'));
par.vartheta=M_.params(strmatch('vartheta',M_.param_names,'exact'));
par.betta=M_.params(strmatch('betta',M_.param_names,'exact'));
par.rho_u=M_.params(strmatch('rho_u',M_.param_names,'exact'));
x_predicted=-par.kappa/(par.kappa^2+par.vartheta*(1-par.betta*par.rho_u));
fprintf('\nAnalytical output response: %10.8f\n',x_predicted)
fprintf('Computed output response: %10.8f\n',oo_.irfs.x_eps_u(1))
figure
subplot(2,2,1)
plot(0:options_.irf-1,oo_.irfs.x_eps_u,'-o')
axis tight
title('Output gap')
legend('Discretion','Location','South')
subplot(2,2,2)
plot(0:options_.irf-1,oo_.irfs.pi_eps_u,'-o')
axis tight
title('Inflation')
subplot(2,2,3)
plot(0:options_.irf-1,oo_.irfs.p_eps_u,'-o')
axis tight
title('Price level')
subplot(2,2,4)
plot(0:options_.irf-1,oo_.irfs.u_eps_u,'-o')
axis tight
title('Cost-push shock')
print -depsc2 Figure_5_2_discretion
