/*
 * This file implements the baseline New Keynesian SOE model of Chapter 8 of
 * Jordi Gali (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton University Press, Second Edition.
 * When the preprocessor variable FDIT is set 1, it replicates Figure 8.1. For each monetary policy rule, it allows
 * replicating the results from Table 8.1. 
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER
 *
 * Notes:
 * - all model variables are expressed in deviations from steady state, i.e. in contrast to
 *      to the chapter, both the nominal interest rate and natural output are not in log-levels, but rather mean 0
 * - In contrast to earlier chapters, eta in Chapter 8 is the trade elasticity, not the semi-elasticity of money demand
 * - In Table 8.1 there seem to be small differences due to rounding
 *
 * This implementation is based on the code by Dmitry Matveev and Johannes Pfeifer. In case you spot mistakes,
 * email us at matveev@uni-mannheim.de or jpfeifer@gmx.de.
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2016-2022 Dmitry Matveev and Johannes Pfeifer
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

//Select monetary policy rule from Chapter 8.5 (FDIT is the baseline Taylor rule used in the chapter) 
@#define FDIT =1

@#define SDIT =0
@#define SCIT =0
@#define PEG  =0
@#define FCIT =0

% define string for saving different policies
@#if FDIT == 1
    case_title='Flexible domestic inflation targeting (FDIT)';    
    short_title='FDIT';
@#else
    @#if SDIT ==1
         case_title='Strict domestic inflation targeting (SDIT)';
        short_title='SDIT';
    @#else
        @#if SCIT ==1
            case_title='Strict CPI inflation targeting (SCIT)';   
            short_title='SCIT';
        @#else
            @#if PEG ==1
                case_title='exchange rate peg (PEG)';
                short_title='PEG';
            @#else
                @#if FCIT ==1
                     case_title='flexible CPI inflation targeting rule (FCIT) (baseline)';
                     short_title='FCIT';
                @#else
                    error('One case must be set to 1')
                @#endif
            @#endif
        @#endif       
    @#endif
@#endif

var y_gap       ${\tilde{y}}$   (long_name='output gap')
    pi_h        ${\pi_H}$       (long_name='domestic inflation')
    i           ${i}$           (long_name='nominal interest rate')
    y_nat       ${y^{n}}$       (long_name='natural output')
    r_nat       ${r^{n}}$       (long_name='natural interest rate')
    s_nat       ${s^{n}}$       (long_name='natural terms of trade')
    y           ${y}$           (long_name='output')
    s_gap       ${\tilde{s}}$   (long_name='terms of trade gap')
    s           ${s}$           (long_name='terms of trade')
    pi          ${\pi}$         (long_name='CPI inflation')
    n           ${n}$           (long_name='employment')
    r_real      ${r^r}$         (long_name='real interest rate')   
    w           ${w}$           (long_name='nominal wage')   
    nx          ${nx}$          (long_name='net exports in terms of domestic outout')   
    c           ${c}$           (long_name='consumption')   
    yhat        ${\hat y}$      (long_name='output deviation from steady state')
    p_h         ${p_H}$         (long_name='domestic price level')
    p           ${p}$           (long_name='CPI')
    er          ${e}$           (long_name='Nominal exchange rate')
    d_er        ${\Delta e}$    (long_name='Nominal exchange rate growth')
    y_star      ${y^*}$         (long_name='world output')
    a           ${a}$           (long_name='AR(1) technology shock process')
    nu          ${\nu}$         (long_name='AR(1) monetary policy shock process')
    z           ${z}$           (long_name='AR(1) preference shock process')
    r_real_ann  ${r^{r,ann}}$   (long_name='annualized real interest rate')
    i_ann       ${i^{ann}}$     (long_name='annualized nominal interest rate')
    r_nat_ann   ${r^{nat,ann}}$ (long_name='annualized natural interest rate')
    pi_ann      ${\pi^{ann}}$   (long_name='annualized CPI inflation rate')
    pi_h_ann    ${\pi_H^{ann}}$ (long_name='annualized domestic inflation rate')
    ;

varexo eps_nu       ${\varepsilon^\nu}$     (long_name='monetary policy shock')
       eps_a        ${\varepsilon^a}$       (long_name='technology shock')
%        eps_y_star   ${\varepsilon^y}$       (long_name='world output growth shock')
       eps_z        ${\varepsilon^z}$       (long_name='preference shock')
       p_star       ${p^*}$                 (long_name='world price level')
       ;

parameters betta    ${\beta}$       (long_name='discount factor')
    siggma          ${\sigma}$      (long_name='inverse EIS')
    varphi          ${\varphi}$     (long_name='inverse Frisch elasticity')
    alppha          ${\alpha}$      (long_name='capital share')
    epsilon         ${\epsilon}$    (long_name='demand elasticity')
    theta           ${\theta}$      (long_name='Calvo parameter')
    upsilon         ${\upsilon}$    (long_name='openness parameter')
    eta             ${\eta}$        (long_name='substitutability foreign/domestic goods')    
    rho_a           ${\rho_a}$      (long_name='autocorrelation technology shock')
    rho_nu          ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')
    rho_y_star      ${\rho_{y^*}}$  (long_name='autocorrelation world output growth shock')
    rho_z           ${\rho_{z}}$    (long_name='autocorrelation preference shock')
    phi_pi          ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    phi_y           ${\phi_{y}}$    (long_name='output feedback Taylor Rule')
    ;

%----------------------------------------------------------------
% Parametrization
%----------------------------------------------------------------
betta   = 0.99;
siggma  = 1;
varphi  = 5;
alppha  = 1/4;
epsilon = 9;
theta   = 3/4;
upsilon = 0.4;
rho_nu  = 0.5;
rho_a   = 0.9;
rho_y_star  = 0;
phi_pi  = 1.5;
phi_y   = 0.5/4;
eta     = 1;
rho_z   =0.5;
%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------
model(linear); 
%Composite parameters
#Omega    =(1-alppha)/(1-alppha+alppha*epsilon);        %p. 233
#lambda   =(1-theta)*(1-betta*theta)/theta*Omega;       %p. 233
#omega    = siggma*eta + (1-upsilon)*(siggma*eta-1);    %p. 235
#Phi = 1 / (1+upsilon*(omega-1));                       %p. 235
#siggma_upsilon = siggma*Phi;                           %p. 235
#Gamma_a  = (1+varphi)/(siggma_upsilon*(1-alppha)+varphi+alppha);                                          %p. 238
#Gamma_star  = - upsilon*(omega-1)*siggma_upsilon*(1-alppha)/(siggma_upsilon*(1-alppha)+varphi+alppha);    %p. 238
#Gamma_z  = - upsilon*omega*Phi*(1-alppha)/(siggma_upsilon*(1-alppha)+varphi+alppha);                      %p. 238   
#kappa_upsilon  = lambda*(siggma_upsilon+(varphi+alppha)/(1-alppha));                                      %p. 238 bottom
#Phi_star=siggma_upsilon*(upsilon*(omega-1)+Gamma_star);                                                %p. 239
#Phi_z=(1-upsilon)*Phi-siggma_upsilon*Gamma_z;                                                          %p. 239        
[name='New Keynesian Phillips Curve (eq. 37)']
pi_h  = betta*pi_h(+1) + kappa_upsilon*y_gap;
[name='Dynamic IS Curve (eq. 29)']
y_gap = y_gap(+1) - 1/siggma_upsilon*(i-pi_h(+1)-r_nat);
[name='Natural output (eq. 35)']
y_nat = Gamma_a*a + Gamma_z*z + Gamma_star*y_star;
[name='Natural rate of interest (eq. 38)']
r_nat = -siggma_upsilon*Gamma_a*(1-rho_a)*a + Phi_star*(y_star(+1)-y_star) + Phi_z*(1-rho_z)*z;
[name='Natural terms of trade (below eq. (35))']
s_nat = siggma_upsilon*(y_nat-y_star)-(1-upsilon)*Phi*z;
[name='Terms of trade gap (middle p. 238)']
s_gap = siggma_upsilon*y_gap;
[name='Output']
y_gap = y - y_nat;
[name=' Terms of trade, p. 238']
s_gap = s - s_nat;
[name='CPI inflation (13)']
pi    = pi_h + upsilon*(s-s(-1));
[name='Production function (eq. 32)']
y     = a + (1-alppha)*n;
[name='Definition real interest rate']
r_real= i - pi_h(+1);
[name='Monetary policy shock, below eq. (39)']
nu    = rho_nu*nu(-1) + eps_nu;
[name='TFP shock, top of p. 233']
a     = rho_a*a(-1) + eps_a;
[name='Preference shock, top of p. 227']
z     = rho_z*z(-1) + eps_z;
[name='FOC wage, eq. (11)']
w-p=siggma*c+varphi*n;
[name='net exports, eq. (31)']
nx=upsilon*(omega/siggma-1)*s-upsilon/siggma*z;
[name='consumption determined by resource constraint, p. 236']
nx=y-c-upsilon*s;
[name='World output growth shock']
%y_star - y_star(-1) = rho_y_star*(y_star(-1) - y_star(-2)) + eps_y_star; %use growth rate rule for world output
y_star  =   0;
[name='Annualized nominal interest rate']
i_ann=4*i;
[name='Annualized real interest rate']
r_real_ann=4*r_real;
[name=' Annualized natural interest rate']
r_nat_ann=4*r_nat;
[name='Annualized CPI inflation']
pi_ann=4*pi;
[name='Annualized domestic inflation']
pi_h_ann=4*pi_h;
[name='Output deviation from steady state']
yhat=y-steady_state(y);
[name='Domestic price level, p. 229']
p_h   = p_h(-1) + pi_h;
[name='CPI definition']
p     = p(-1) + pi;
[name='Nominal exchange rate']
s     = er + p_star -  p_h ;
[name='Definiion exchange rate growth']
d_er=er-er(-1);

//Set policy rule according to Chapter 8.5                
@#if FDIT == 1
    [name='Interest Rate Rule (eq. 39) (FDIT)']
    i  = phi_pi*pi_h + phi_y*yhat + nu;
@#else
    @#if SDIT ==1
        [name='SDIT']
        pi_h=0;
    @#else
        @#if SCIT ==1
            [name='SCIT']
            pi=0;   
        @#else
            @#if PEG ==1
                [name='PEG']
                er=0;
            @#else
                @#if FCIT ==1
                    [name='FCIT']
                    i  = phi_pi*pi + phi_y*yhat + nu;    
                @#endif
            @#endif
        @#endif       
    @#endif
@#endif
end;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------
shocks;
    var eps_nu = 0.25^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
end;

%----------------------------------------------------------------
%  steady states: all 0 due to linear model
%---------------------------------------------------------------
steady;
check;

%----------------------------------------------------------------
% generate IRFs, replicates Figures 8.1, p. 243
%----------------------------------------------------------------
stoch_simul(order = 1,irf=15) y_gap n pi_h_ann pi_ann p_h p s er i_ann r_real_ann;


%----------------------------------------------------------------
% replicate Table 8.1, p. 251
%----------------------------------------------------------------

%set parameters to the special case described
set_param_value('siggma',1);
set_param_value('eta',1);
set_param_value('upsilon',0.4);
shocks;
    var eps_a = 1; //1 percent shock
    var eps_nu = 0; 
%     var eps_y_star = 0;
    var eps_z = 0;
    var p_star = 0;
end;

%compute variances for relevant variables
stoch_simul(order = 1,irf=0) y y_gap pi_h pi s d_er ;
%find location of variances in results
y_pos=strmatch('y',var_list_ ,'exact');
y_gap_pos=strmatch('y_gap',var_list_ ,'exact');
pi_pos=strmatch('pi',var_list_ ,'exact');
pi_h_pos=strmatch('pi_h',var_list_ ,'exact');
s_pos=strmatch('s',var_list_ ,'exact');
d_er_pos=strmatch('d_er',var_list_ ,'exact');

%get recent parameters for computation
par.upsilon=M_.params(strmatch('upsilon',M_.param_names,'exact'));
par.theta=M_.params(strmatch('theta',M_.param_names,'exact'));
par.alppha=M_.params(strmatch('alppha',M_.param_names,'exact'));
par.betta=M_.params(strmatch('betta',M_.param_names,'exact'));
par.epsilon=M_.params(strmatch('epsilon',M_.param_names,'exact'));
par.siggma=M_.params(strmatch('siggma',M_.param_names,'exact'));
par.varphi=M_.params(strmatch('varphi',M_.param_names,'exact'));
par.lambda=(1-par.theta)*(1-par.betta*par.theta)/par.theta*(1-par.alppha)/(1-par.alppha+par.alppha*par.epsilon);

%read out variances
variance.y=oo_.var(y_pos,y_pos);
variance.y_gap=oo_.var(y_gap_pos,y_gap_pos);
variance.pi_h=oo_.var(pi_h_pos,pi_h_pos);
variance.pi=oo_.var(pi_pos,pi_pos);
variance.s=oo_.var(s_pos,s_pos);
variance.d_er=oo_.var(d_er_pos,d_er_pos);
%compute loss function
L=(1-par.upsilon)*0.5*(((1+par.varphi)/(1-par.alppha))*variance.y_gap+par.epsilon/par.lambda*variance.pi_h)/100;

//Print result
labels={'sigma(y)';'sigma(tilde y)';'sigma(pi_h)';'sigma(pi)';'sigma(s)';'sigma(Delta e)';'L'};
headers={' ',short_title};
values=[sqrt([variance.y;variance.y_gap;variance.pi_h;variance.pi;variance.s;variance.d_er]);L];
options_.noprint=0;
dyntable(options_,'Table 8.1: Properties of simple policy rules',headers,labels,values,size(labels,2)+2,4,3)

