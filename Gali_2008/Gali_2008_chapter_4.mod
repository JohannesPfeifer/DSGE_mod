/*
 * This file performs welfare analysis of Chapter 4.4 on simple rules in the baseline New Keynesian model
 * of Jordi Gali (2008): Monetary Policy, Inflation, and the Business Cycle, Princeton University Press.
 * It replicates Table 4.1 (Taylor Rule and Constant Money Growth) .
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER
 *
 * Notes:
 * - all model variables are expressed in deviations from steady state, i.e. in contrast to
 *      to the chapter, both the nominal interest rate and natural output are not in log-levels, but rather mean 0
 *
 * This implementation is based on the code by Dmitry Matveev and Johannes Pfeifer. In case you spot mistakes,
 * email us jpfeifer@gmx.de.
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2016 Dmitry Matveev and Johannes Pfeifer
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

/*
 * Derived from Gali_2015_chapter_4.mod
 * modified by Detlef Hering to reproduce table 4.1 (page 83) of Gali 2008
 */

%define whether to use interest rate or money growth rate rule 
@#ifndef money_growth_rule
@#define money_growth_rule=1
@#endif

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
    @#if money_growth_rule==0
        nu      ${\nu}$                 (long_name='AR(1) monetary policy shock process')    
    @#else
        money_growth  ${\Delta m_q}$      (long_name='money growth')
        zeta          ${\zeta}$           (long_name='level of money demand shock')
		Delta_zeta    ${\Delta \zeta}$    (long_name='Difference of money demand shock')
   @#endif
    a           ${a}$                   (long_name='AR(1) technology shock process')
    r_real_ann  ${r^{r,ann}}$           (long_name='annualized real interest rate')
    i_ann       ${i^{ann}}$             (long_name='annualized nominal interest rate')
    r_nat_ann   ${r^{nat,ann}}$         (long_name='annualized natural interest rate')
    pi_ann      ${\pi^{ann}}$           (long_name='annualized inflation rate')
    z           ${z}$                   (long_name='AR(1) preference shock process')
;     

varexo  eps_a       ${\varepsilon_a}$       (long_name='technology shock')
        @#if money_growth_rule==0
            eps_nu  ${\varepsilon_\nu}$     (long_name='monetary policy shock')
        @#else   
            eps_zeta       ${\varepsilon_\zeta}$   (long_name='money demand shock innovation')
        @#endif
        eps_z       ${\varepsilon_z}$   (long_name='preference shock innovation')
       ;

parameters alppha       ${\alppha}$     (long_name='capital share')
    betta               ${\beta}$       (long_name='discount factor')
    rho_a               ${\rho_a}$      (long_name='autocorrelation technology shock')
    @#if money_growth_rule==0
        rho_nu          ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')
    @#else   
        rho_zeta        ${\rho_{\zeta}}$ (long_name='autocorrelation monetary demand shock')
		sigma_zeta      ${\sigma_{\zeta}}$ (long_name='inverse EIS monetary demand shock')
    @#endif
    rho_z               ${\rho_{z}}$    (long_name='autocorrelation monetary demand shock')
    siggma              ${\sigma}$      (long_name='inverse EIS')
    varphi              ${\varphi}$     (long_name='inverse Frisch elasticity')
    phi_pi              ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    phi_y               ${\phi_{y}}$    (long_name='output feedback Taylor Rule')
    eta                 ${\eta}$        (long_name='semi-elasticity of money demand')
    epsilon             ${\epsilon}$    (long_name='demand elasticity')
    theta               ${\theta}$      (long_name='Calvo parameter')
	rho                 ${\rho}$        (long_name='interest rate parameter')
    ;
%----------------------------------------------------------------
% Parametrization, p. 52, 55
%----------------------------------------------------------------
siggma = 1;
varphi=1;
phi_pi = 1.5;
phi_y  = 0.125;
theta=2/3;
@#if money_growth_rule==0
    rho_nu =0.5;
@#else   
    rho_zeta=0.6;
	sigma_zeta=0.63;
@#endif
rho_z=0.5;
rho_a  = 0.9;
betta = 0.99;
eta  = 4;
alppha=1/3;
epsilon=6;
rho = -log(betta);

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 
//Composite parameters
#Omega=(1-alppha)/(1-alppha+alppha*epsilon);        %defined on page 47
#psi_n_ya=(1+varphi)/(siggma*(1-alppha)+varphi+alppha);   %defined on page 48
#lambda=(1-theta)*(1-betta*theta)/theta*Omega;      %defined on page 47
#kappa=lambda*(siggma+(varphi+alppha)/(1-alppha));     %defined on page 49
[name='New Keynesian Phillips Curve eq. (22)']
pi=betta*pi(+1)+kappa*y_gap;
[name='Dynamic IS Curve eq. (23)']
y_gap=-1/siggma*(i-pi(+1)-r_nat)+y_gap(+1);
@#if money_growth_rule==0
    [name='Interest Rate Rule eq. (26)']
    i=rho+phi_pi*pi+phi_y*yhat+nu;
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
    [name='8. Monetary policy shock']
    nu=rho_nu*nu(-1)+eps_nu;
@#endif
[name='TFP shock']
a=rho_a*a(-1)+eps_a;
[name='Production function (eq. 14)']
y=a+(1-alppha)*n;
[name='Preference shock, p. 54']
z     = rho_z*z(-1) + eps_z;
[name='Money growth (derived from eq. (4))']
m_growth_ann=4*(y-y(-1)-eta*(i-i(-1))+pi);
@#if money_growth_rule==1
    [name='Real money demand, (p. 114)']
    m_real=y-eta*i - zeta;
    [name='definition nominal money growth']
    money_growth=m_real-m_real(-1)+pi;
    [name='0 exogenous money growth (p. 114)']
    money_growth=0;
    [name='Definition of first difference of money demand shock']
	Delta_zeta = zeta -zeta(-1);
    [name='exogenous process for money demand shock']
	Delta_zeta = rho_zeta*Delta_zeta(-1)+sigma_zeta*eps_zeta;
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

end;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------

@#if money_growth_rule==0
    shocks;
        var eps_a  = 1^2; //unit shock to technology
    end;
@#else
    shocks;
        var eps_a  = 1^2; //unit shock to technology
        var eps_zeta  = 1^2;        
    end;
@#endif

   

%----------------------------------------------------------------
%  steady states: all 0 due to linear model in log deviations
%---------------------------------------------------------------
%
steady;
check;

%----------------------------------------------------------------
% solve for welfare-relevant second moments (technology shocks) 
% replicates part of Table 4.1, p. 113 (interest rate rule)
%                          4.2, p. 119 (money growth rule)
%----------------------------------------------------------------
@#if money_growth_rule==0
	stoch_simul(order=1,irf=0,noprint) y y_gap pi pi_ann a;
@#else
	stoch_simul(order=1,irf=0,noprint) y y_gap pi pi_ann m_real a z zeta Delta_zeta;
@#endif

%----------------------------------------------------------------
% get variable positions in variable list
%----------------------------------------------------------------

y_gap_pos=strmatch('y_gap',var_list_ ,'exact');
pi_pos=strmatch('pi',var_list_ ,'exact');


@#if money_growth_rule==0

    %----------------------------------------------------------------
    % Replicate panel a) of Table 4.1 (Taylor Rule)
    %----------------------------------------------------------------

    phi_pi_vec=[1.5 1.5 5 1.5];
    phi_y_vec=[0.125 0 0 1];
	
    variance.y_gap=NaN(1,length(phi_pi_vec));
    variance.y=NaN(1,length(phi_pi_vec));
    variance.pi=NaN(1,length(phi_pi_vec));
    L=NaN(1,length(phi_pi_vec));
    for ii=1:length(phi_pi_vec)
        set_param_value('phi_pi',phi_pi_vec(ii));
        set_param_value('phi_y',phi_y_vec(ii));
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %loop over stoch_simul
        if ~info(1)
            %read out current parameter values
            par.theta=M_.params(strmatch('theta',M_.param_names,'exact'));
            par.alppha=M_.params(strmatch('alppha',M_.param_names,'exact'));
            par.betta=M_.params(strmatch('betta',M_.param_names,'exact'));
            par.epsilon=M_.params(strmatch('epsilon',M_.param_names,'exact'));
            par.siggma=M_.params(strmatch('siggma',M_.param_names,'exact'));
            par.varphi=M_.params(strmatch('varphi',M_.param_names,'exact'));
            
            par.lambda=(1-par.theta)*(1-par.betta*par.theta)/par.theta*(1-par.alppha)/(1-par.alppha+par.alppha*par.epsilon);
            
            variance.y_gap(ii)=oo_.var(y_gap_pos,y_gap_pos);
            variance.pi(ii)=oo_.var(pi_pos,pi_pos);
            L(ii)=0.5*((par.siggma+(par.varphi+par.alppha)/(1-par.alppha))*variance.y_gap(ii)+par.epsilon/par.lambda*variance.pi(ii))/100;
        end
    end
    //Print result

    labels={'phi_pi';'phi_y';'sigma(tilde y)';'sigma(pi)';'L'};
    headers={' ';' ';' ';' '};	
    values=[phi_pi_vec;phi_y_vec;sqrt(variance.y_gap);4*sqrt(variance.pi);L];
    options_.noprint=0;
    dyntable(options_,'Taylor Rule',headers,labels,values,size(labels,2)+2,4,3)
    options_.noprint=1;

@#else
    rho_zeta_vec=[0 0.6];
    sigma_zeta_vec=[0 0.63];
    variance.y_gap=NaN(1,length(rho_zeta_vec));
    variance.y=NaN(1,length(rho_zeta_vec));
    variance.pi=NaN(1,length(rho_zeta_vec));
    L=NaN(1,length(rho_zeta_vec));
    for ii=1:length(rho_zeta_vec)
        set_param_value('rho_zeta',rho_zeta_vec(ii));
        set_param_value('sigma_zeta',sigma_zeta_vec(ii));
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %loop over stoch_simul
        if ~info(1)
            %read out current parameter values
            par.theta=M_.params(strmatch('theta',M_.param_names,'exact'));
            par.alppha=M_.params(strmatch('alppha',M_.param_names,'exact'));
            par.betta=M_.params(strmatch('betta',M_.param_names,'exact'));
            par.epsilon=M_.params(strmatch('epsilon',M_.param_names,'exact'));
            par.siggma=M_.params(strmatch('siggma',M_.param_names,'exact'));
            par.varphi=M_.params(strmatch('varphi',M_.param_names,'exact'));
            
            par.lambda=(1-par.theta)*(1-par.betta*par.theta)/par.theta*(1-par.alppha)/(1-par.alppha+par.alppha*par.epsilon);
            
            variance.y_gap(ii)=oo_.var(y_gap_pos,y_gap_pos);
            variance.pi(ii)=oo_.var(pi_pos,pi_pos);
            L(ii)=0.5*((par.siggma+(par.varphi+par.alppha)/(1-par.alppha))*variance.y_gap(ii)+par.epsilon/par.lambda*variance.pi(ii))/100;
        end
    end
    //Print result

    labels={'rho_zeta';'sigma_zeta';'sigma(tilde y)';'sigma(pi)';'L'};
    headers={' ';' ';' ';' '};	
    values=[rho_zeta_vec;sigma_zeta_vec;sqrt(variance.y_gap);4*sqrt(variance.pi);L];
    options_.noprint=0;
    dyntable(options_,'Constant money growth',headers,labels,values,size(labels,2)+2,4,3)
@#endif
