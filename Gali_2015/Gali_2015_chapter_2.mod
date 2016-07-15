/*
 * This file implements the baseline Classical Monetary Economy model of Jordi Galí (2015): Monetary Policy, Inflation,
 * and the Business Cycle, Princeton University Press, Second Edition, Chapter 2
 *
 * Note that this mod-file implements the non-linear first order conditions and that the IRFs show the linear deviations
 * from steady state.
 *
 * It demonstrate the neutrality of money by showing that real variables do not move after a monetary policy shock
 * 
 * THIS MOD-FILE REQUIRES DYNARE 4.5 (I.E. THE CURRENT UNSTABLE VERSION)
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

var C               ${C}$           (long_name='Consumption')
    W_real          ${\frac{W}{P}}$ (long_name='Real Wage')
    Pi              ${\Pi}$         (long_name='inflation')
    A               ${A}$           (long_name='AR(1) technology process')
    N               ${N}$           (long_name='Hours worked')
    R               ${R^n}$         (long_name='Nominal Interest Rate') 
    realinterest    ${R^{r}}$       (long_name='Real Interest Rate')
    Y               ${Y}$           (long_name='Output') 
    nu              ${\nu}$         (long_name='AR(1) monetary policy shock process')    
    m_growth_ann    ${\Delta M}$    (long_name='money growth')
    Q               ${Q}$           (long_name='Bond price')
    Z               ${Z}$           (long_name='AR(1) preference shock process')
    ;

varexo eps_a    ${\varepsilon_a}$   (long_name='technology shock')
       eps_z    ${\varepsilon_z}$   (long_name='preference shock')
       eps_nu   ${\varepsilon_\nu}$ (long_name='monetary policy shock')
       ;   

parameters alppha   ${\alpha}$      (long_name='capital share')
    betta           ${\beta}$       (long_name='discount factor')
    rho_a           ${\rho_a}$      (long_name='autocorrelation technology shock')
    rho_z           ${\rho_z}$      (long_name='autocorrelation preference shock')
    rho_nu          ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')
    siggma          ${\sigma}$      (long_name='log utility')
    varphi          ${\varphi}$     (long_name='unitary Frisch elasticity')
    phi_pi          ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    eta             ${\eta}$        (long_name='semi-elasticity of money demand')
    ;

%----------------------------------------------------------------
% Follows parametrization of Chapter 3, p. 52
%----------------------------------------------------------------

alppha = 1/4;
betta  = 0.99;
rho_a  = 0.9;
rho_z  = 0.5;
rho_nu = 0.5;
siggma = 1;
varphi = 5;
phi_pi = 1.5;
eta    = 3.77; %footnote 11, p. 115


%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model;
[name='FOC Wages, eq. (7)']
W_real=C^siggma*N^varphi;
[name='Euler equation eq. (8)']
Q=betta*(C(+1)/C)^(-siggma)*(Z(+1)/Z);
[name='Definition nominal interest rate), p. 22 top']
R=1/Q;
[name='Production function eq. (12)']
Y=A*N^(1-alppha);
[name='FOC wages firm, eq. (14)']
W_real=(1-alppha)*A*N^(-alppha);
[name='Definition Real interest rate, eq. 22']
R=realinterest*Pi(+1);
[name='Monetary Policy Rule, p. 26 bottom/eq. (22)']
R=1/betta*Pi^phi_pi*exp(nu);
[name='Market Clearing, eq. (15)']
C=Y;
[name='Technology Shock, p.22']
log(A)=rho_a*log(A(-1))+eps_a;
[name='Preference Shock, p.21']
log(Z)=rho_z*log(Z(-1))+eps_z;
[name='Monetary policy shock, p.23']
nu=rho_nu*nu(-1)+eps_nu;
[name='Money growth (derived from eq. (11))']
m_growth_ann=4*(log(C)-log(C(-1))-eta*(log(R)-log(R(-1)))+log(Pi));
end;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------

shocks;
var eps_a; stderr 1;
var eps_z; stderr 1;
var eps_nu; stderr 1;
end;

%----------------------------------------------------------------
%  Initial Values for steady state
%---------------------------------------------------------------

steady_state_model;
A=1;
Z=1;
R=1/betta;
Pi=1;
Q=1/R;
realinterest=R;
N=(1-alppha)^(1/((1-siggma)*alppha+varphi+siggma));
C=A*N^(1-alppha);
W_real=(1-alppha)*A*N^(-alppha);
Y=C;
m_growth_ann=0;
end;


resid(1);
steady;
check;

%----------------------------------------------------------------
% generate IRFs to show neutrality of money
%----------------------------------------------------------------
write_latex_dynamic_model;
stoch_simul(irf=20,order=1) Y C Pi R realinterest m_growth_ann;