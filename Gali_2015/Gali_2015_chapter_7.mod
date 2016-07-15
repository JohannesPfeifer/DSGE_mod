/*
 * This file implements the New Keynesian model with price and wage rigidities and unemployment 
 * of Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton 
 * University Press, Second Edition, Chapter 7
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5 (I.E. THE CURRENT UNSTABLE VERSION)
 *
 * Notes:
 *  - all model variables are expressed in deviations from steady state, i.e. in contrast to
 *      to the chapter, the nominal interest rate, natural output, the natural real wage, the unemployment rate,
 *      and the labor force are not in log-levels, but rather mean 0
 *  - on page 212 it is noted that "the optimality conditions are identical to those derived in Chapter 3". However,
 *      it should be "identical to Chapter 6" as Chapter 3 does not feature wage stickiness. This can be verified 
 *      by comparing Chapter 7, equations (18)-(21), p. 212 to the identical ones in Chapter 6, eq. (30)-(33) on p. 181
 *  - To replicate Figures 7.5 and 7.6, run all three policy cases by setting the respective indicators 
 *      (taylor_rule,optimal_policy,simple_rule_unemployment) in turn to 1. When the code detects all three 
 *      saved results files, the figure will be displayed   
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
@#define taylor_rule=0
@#define optimal_policy=1
@#define simple_rule_unemployment=0
        
@#if taylor_rule
    save_name='IRFs_taylor_rule';
@#else
    @#if simple_rule_unemployment
        save_name='IRFs_simple_rule_unemployment';
    @#else
        @#if optimal_policy
            save_name='IRFs_optimal_policy';
        @#else
            error('One dummy must be set to 1')
        @#endif 
    @#endif 
@#endif 
            
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
    u           ${u}$                   (long_name='unemployment')
    l           ${l}$                   (long_name='labor force')
    uhat        ${\hat u}$              (long_name='unemployment gap')
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
    u_n                 ${u^n}$         (long_name='natural rate of unemployment')
    ;
%----------------------------------------------------------------
% Parametrization, p. 67, p. 113-115, and p. 208
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
u_n=0;  %normalized to 0 here as all variables are demeaned
epsilon_w=4.5;
theta_w=3/4;
%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 
//Composite parameters
#Omega=(1-alppha)/(1-alppha+alppha*epsilon_p);              %defined on page 166
#psi_n_ya=(1+varphi)/(siggma*(1-alppha)+varphi+alppha);     %defined on page 208
#psi_n_wa=(1-alppha*psi_n_ya)/(1-alppha);                   %defined on page 208
#aleph_p=alppha*lambda_p/(1-alppha);                         %defined on page 172
#aleph_w=lambda_w*(siggma+varphi/(1-alppha));               %defined on page 172
[name='New Keynesian Phillips Curve eq. (13)']
pi_p=betta*pi_p(+1)+aleph_p*y_gap+lambda_p*w_gap;
[name='New Keynesian Wage Phillips Curve eq. (14) (see footnote 8, p. 208)']
pi_w=betta*pi_w(+1)+aleph_w*y_gap-lambda_w*w_gap;
[name='Dynamic IS Curve eq. (12)']
y_gap=-1/siggma*(i-pi_p(+1)-r_nat)+y_gap(+1);
@#if taylor_rule
    [name='Interest Rate Rule eq. (17)']
    i=phi_pi*pi_p+phi_y*yhat+nu;
@#else
    @#if simple_rule_unemployment
        [name='Interest Rate Rule eq. (17)']
        i=1.5*pi_p-0.5*uhat+nu;
    @#else
        %ramsey policy, no rule required
    @#endif 
@#endif 
        

[name='Definition natural rate of interest, p. 208']
r_nat=-siggma*psi_n_ya*(1-rho_a)*a+(1-rho_z)*z;
[name='Definition wage gap, Ch. 6, eq (21)']
w_gap=w_gap(-1)+pi_w-pi_p-(w_nat-w_nat(-1));
[name='Definition natural wage, p. 208']
w_nat=psi_n_wa*a;
[name='Definition natural wage, p. 208']
varphi*uhat=w_gap-(siggma+varphi/(1-alppha))*y_gap;
[name='definition unemployment, eq. (7)']
uhat=u-u_n;
u=l-n;
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
z     = rho_z*z(-1) + eps_z;
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


@#if taylor_rule

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------

shocks;
    var eps_nu = 0.25^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
end;

%----------------------------------------------------------------
% generate IRFs for monetary policy shock, replicates Figures 7.2, p. 209
%----------------------------------------------------------------
stoch_simul(order = 1,irf=15) u n l w_real;

stoch_simul(order = 1,irf=0,noprint,contemporaneous_correlation) u y;

%----------------------------------------------------------------
% replicate Table 7.1, p. 211
%----------------------------------------------------------------

theta_w_vec=[0.1 0.5 0.75];
rho_nu_vec=[0 0.5 0.9];
u_pos=strmatch('u',var_list_,'exact');
y_pos=strmatch('y',var_list_,'exact');

verbatim;
volatility_mat=NaN(length(rho_nu_vec),length(theta_w_vec));
correlation_mat=NaN(length(rho_nu_vec),length(theta_w_vec));
cyclicality_mat=NaN(length(rho_nu_vec),length(theta_w_vec));

for theta_w_iter=1:length(theta_w_vec)
    for rho_nu_iter=1:length(rho_nu_vec)
        set_param_value('theta_w',theta_w_vec(theta_w_iter));
        set_param_value('rho_nu',rho_nu_vec(rho_nu_iter));
        info=stoch_simul(var_list_);
        volatility_mat(rho_nu_iter,theta_w_iter)=sqrt(oo_.var(u_pos,u_pos));
        correlation_mat(rho_nu_iter,theta_w_iter)=oo_.autocorr{1,1}(u_pos,u_pos);
        cyclicality_mat(rho_nu_iter,theta_w_iter)=oo_.contemporaneous_correlation(u_pos,y_pos);
    end
end
options_.noprint=0;
labels=strvcat('sigma(y)','sigma(tilde y)','sigma(pi)','L');
headers=strvcat('theta_v:', num2str(theta_w_vec'));
labels=[repmat('rho_nu=',length(rho_nu_vec),1),num2str(rho_nu_vec')];
dyntable(options_,'Volatility',headers,labels,volatility_mat,size(labels,2)+2,4,3)
    
dyntable(options_,'Persistence',headers,labels,correlation_mat,size(labels,2)+2,4,3)

dyntable(options_,'Cyclicality',headers,labels,cyclicality_mat,size(labels,2)+2,4,3)

end;
@#endif

        
%reset parameters
set_param_value('theta_w',3/4);
set_param_value('rho_nu',0.5);

shocks;
    var eps_nu = 0; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    var eps_a = 1;
    var eps_z = 1;
end;

@#if optimal_policy
//planner objective, uses lambda_w and lambda_p updated in steady_state_model-block
    planner_objective 0.5*((siggma+(varphi+alppha)/(1-alppha))*y_gap^2+ epsilon_p/lambda_p*pi_p^2+epsilon_w*(1-alppha)/lambda_w*pi_w^2);
    ramsey_policy(instruments=(i),irf=17,planner_discount=betta,noprint) y u n l w_real pi_p_ann;
@#else    
    stoch_simul(order=1,irf=17) y u n l w_real pi_p_ann;
@#endif
save(save_name,'oo_','var_list_');


%----------------------------------------------------------------
% replicate Figures 7.4 and 7.5, p. 215-16 (requires the file to 
%   be run with all three policy choices first)
%----------------------------------------------------------------

if exist('IRFs_taylor_rule.mat','file') && exist('IRFs_simple_rule_unemployment.mat','file') && exist('IRFs_optimal_policy.mat')
    IRFs_taylor_rule=load('IRFs_taylor_rule');
    IRFs_simple_rule_unemployment=load('IRFs_simple_rule_unemployment');
    IRFs_optimal_policy=load('IRFs_optimal_policy');
    figure('Name','Optimal policy vs simple rule: technology shock')
    iter=1;
    shock_name='a';
    IRF_length=length(IRFs_taylor_rule.oo_.irfs.([deblank(var_list_(1,:)),'_eps_',shock_name]));
    for var_name_iter=1:size(var_list_,1)
        subplot(3,2,iter)
        plot(1:IRF_length,IRFs_taylor_rule.oo_.irfs.([deblank(var_list_(var_name_iter,:)),'_eps_',shock_name]),'-s',1:IRF_length,IRFs_simple_rule_unemployment.oo_.irfs.([deblank(var_list_(var_name_iter,:)),'_eps_',shock_name]),'-o',1:IRF_length,IRFs_optimal_policy.oo_.irfs.([deblank(var_list_(var_name_iter,:)),'_eps_',shock_name]),'-d')
        axis tight
        iter=iter+1;
    end
    figure('Name','Optimal policy vs simple rule: demand shock')
    iter=1;
    shock_name='z';
    IRF_length=length(IRFs_taylor_rule.oo_.irfs.([deblank(var_list_(1,:)),'_eps_',shock_name]));
    for var_name_iter=1:size(var_list_,1)
        subplot(3,2,iter)
        plot(1:IRF_length,IRFs_taylor_rule.oo_.irfs.([deblank(var_list_(var_name_iter,:)),'_eps_',shock_name]),'-s',1:IRF_length,IRFs_simple_rule_unemployment.oo_.irfs.([deblank(var_list_(var_name_iter,:)),'_eps_',shock_name]),'-o',1:IRF_length,IRFs_optimal_policy.oo_.irfs.([deblank(var_list_(var_name_iter,:)),'_eps_',shock_name]),'-d')
        axis tight
        iter=iter+1;
    end
end