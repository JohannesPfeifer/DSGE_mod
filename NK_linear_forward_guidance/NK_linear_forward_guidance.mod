/*
 * This file implements forward guidance in the the baseline New Keynesian model of Jordi Galí (2015): Monetary Policy, Inflation,
 * and the Business Cycle, Princeton University Press, Second Edition, Chapter 3
 *
 * The model employs a perfect foresight setup and computes the required values of a monetary policy shock 
 * to achieve a given path of the nominal interest rate: fixed at the steady state for three periods, followed 
 * by a drop of 100 basis points annualized (25bp quarterly).
 *
 * THIS MOD-FILE REQUIRES DYNARE 6.0 OR HIGHER
 *
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2021 Johannes Pfeifer
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
    nu          ${\nu}$                 (long_name='AR(1) monetary policy shock process')        
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
        eps_nu      ${\varepsilon_\nu}$     (long_name='monetary policy shock')
        eps_z       ${\varepsilon_z}$   (long_name='preference shock innovation')
       ;

parameters alppha       ${\alpha}$     (long_name='capital share')
    betta               ${\beta}$       (long_name='discount factor')
    rho_a               ${\rho_a}$      (long_name='autocorrelation technology shock')
    rho_nu              ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')   
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
% Parametrization, p. 67  and p. 113-115
%----------------------------------------------------------------
% Parametrization, p. 67-75
%----------------------------------------------------------------
siggma = 1;         %p. 67
varphi=5;           %p. 67, Frisch elasticity of 0.2
phi_pi = 1.5;       %p. 68 
phi_y  = 0.125;     %p. 68 (5/4)
theta=3/4;          %p. 67
rho_nu =0.5;    %p. 68
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
[name='Interest Rate Rule eq. (26)']
i=phi_pi*pi+phi_y*yhat+nu;
[name='Definition natural rate of interest eq. (24)']
r_nat=-siggma*psi_n_ya*(1-rho_a)*a+(1-rho_z)*z;
[name='Definition real interest rate']
r_real=i-pi(+1);
[name='Definition natural output, eq. (20)']
y_nat=psi_n_ya*a;
[name='Definition output gap']
y_gap=y-y_nat;
[name='Monetary policy shock']
nu=rho_nu*nu(-1)+eps_nu;
[name='TFP shock']
a=rho_a*a(-1)+eps_a;
[name='Production function (eq. 14)']
y=a+(1-alppha)*n;
[name='Preference shock, p. 54']
z     = rho_z*z(-1) - eps_z;
[name='Money growth (derived from eq. (4))']
m_growth_ann=4*(y-y(-1)-eta*(i-i(-1))+pi);
[name='Real money demand (eq. 4)']
m_real=y-eta*i;
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
%  steady states: all 0 due to linear model
%---------------------------------------------------------------
resid;
steady;
check;

perfect_foresight_setup(periods=100);

options_.verbosity=0;
start_value=[0,0,0,0]';
shock_name='eps_nu';
target_value=[0,0,0,-1];
target_name='i_ann';

% %available in Dynare 5.0
% options_.jacobian_flag=0;
% [x, errorflag, fvec, fjac] = dynare_solve('distance',start_value, options_, shock_name,target_value,target_name,M_,options_,oo_);
[x, errorflag] = csolve('distance',start_value,[],1e-6,500,shock_name,target_value,target_name,M_,options_,oo_);


shocks;
var eps_nu;
periods 1:4;
values (x);
end;


perfect_foresight_setup(periods=100);
perfect_foresight_solver;


figure('Name','Time Path')
T=10;
subplot(3,1,1)
plot(oo_.endo_simul(strmatch('y',M_.endo_names,'exact'),1:M_.maximum_lag+T),'-x')
title('Output')
ylabel('% dev. from SS')

subplot(3,1,2)
plot(100*oo_.endo_simul(strmatch('i_ann',M_.endo_names,'exact'),1:M_.maximum_lag+T),'-x')
title('Annualized Nominal Interest Rate')
ylabel('Dev. from SS (basis points)')

subplot(3,1,3)
plot(100*oo_.exo_simul(1:M_.maximum_lag+T,strmatch('eps_nu',M_.exo_names,'exact')),'-x')
title('Monetary policy shock')
ylabel('Basis points')

