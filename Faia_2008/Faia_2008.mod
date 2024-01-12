/*
 * This files tries and fails to replicate the results of Faia (2008): Optimal monetary policy rules with
 *     labor market frictions, Journal of Economic Dynamics & Control 32 (2008) 1600–1621
 * 
 *  Notes:
 *      - The IRFs in Figures 1 and 2 for the Tp case (presumably the strong inflation response) are not obtained 
 *          using the present codes. It's unclear what causes this discrepancy. A few potential causes are 
 *          outlined below:
 *      - The paper fails to clearly distinguish between gross output produced by firms 
 *         before subtracting adjustment and vacany costs and net output after subtracting those
 *         For example, equation (15) should depend on gross output while equations (26)-(27) should 
 *         depend on net output. It's not clear what was used in the original codes.
 *      - The stated 6% unemployment rate (p. 1609: the steady-state unemployment rate is set equal to 0.06) cannot 
 *          be a calibration target. Rather, its value is alread implied by the other targets for q, theta*q, and rho.
 *          The present file obtains a steady state unemployment rate of 12.7%, given those targets.
 *      - It is not clear what the units in Figures 1 and 2 are and what the shock size was. I assume that all IRFs show
 *          percentage deviations from steady state (i.e. variables are in logs) and that a unit shock was used. That roughly 
 *          matches the qualitative and quantative results.
 *      - Equation (16): It must be mu_t - (kappa/q) instead of mu_{t-}*(kappa/q)
 *      - The AR process below equation (6) should have *exp(epsilon_t) 
*/

/*
 * Copyright (C) 2023 Johannes Pfeifer
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
@#ifndef Ramsey
    @#define Ramsey=0
@#endif


var lambda  ${\Lambda}$     (long_name='Lagrange multiplier')
    c       ${c}$           (long_name='Consumption')
    R       ${R}$           (long_name='Nominal interest rate')
    pi      ${\pi}$         (long_name='Inflation rate')
    theta   ${\theta}$      (long_name='market tightness')
    v       ${v}$           (long_name='vacancies')
    u       ${u}$           (long_name='Unemployment rate')
    m       ${m}$           (long_name='matches')
    q       ${q}$           (long_name='Meeting rate between firms and workers') 
    n       ${n}$           (long_name='employment')
    y_gross ${y^{gross}}$   (long_name='gross output')
    y_net   ${y^{net}}$     (long_name='gross output')
    mu      ${\mu}$         (long_name='Lagrange multiplier')
    z       ${z}$           (long_name='log TFP')
    mc      ${mc}$          (long_name='marginal costs')
    w       ${w}$           (long_name='real wage')
    g       ${g}$           (long_name='government spending')
    z_G     ${g}$           (long_name='government spending shock')
    log_y_net   ${\log y}$      (long_name='log output')
    log_v       ${\log v}$      (long_name='log vacancies')
    log_w       ${\log w}$      (long_name='log wages')
    log_u       ${\log u}$      (long_name='log unemployment')
    log_theta   ${\log \theta}$      (long_name='log tightness')
    log_pi      ${\log \pi}$    (long_name='log tightness')
    ;

varexo epsilon_G    ${\varepsilon_G}$   (long_name='Government spending shock')
    epsilon_z       ${\varepsilon_Z}$   (long_name='Technology shock')
    ;

PARAMETERS 
    epsilon     ${\varepsilon}$         (long_name='substitution elasticity')
    Psi         ${\psi}$                (long_name='Price adjustment costs')      
    betta       ${\beta}$               (long_name='discount factor')      
    xi          ${\xi}$                 (long_name='exponent matching function')     
    varsigma    ${\varsigma}$           (long_name='bargaining power') 
    rho         ${\rho}$                (long_name='separation rate') 
    m_param     ${m^{par}}$             (long_name='scaling parameter in matching function (5)')
    b_w_target  ${\frac{b}{\bar w}}$    (long_name='target value of steady state b/w')
    b           $b$                     (long_name='real unemployment benefits') 
    kappa       ${\kappa}$              (long_name='vacancy posting cost') 
    lambda_par  ${\lambda}$             (long_name='wage rigidity')
    g_share     ${\frac{\bar G}{\bar Y}}$          (long_name='steady state government spending share') 
    G_SS        ${\bar G}$              (long_name='steady state government spending') 
    rho_G       ${\rho_G}$              (long_name='persistence government spending')
    rho_Z       ${\rho_Z}$              (long_name='persistence TFP')
    siggma      ${\sigma}$              (long_name='risk aversion')
    @#if Ramsey==0
    phi_R       ${\phi_R}$              (long_name='interest rate smoothing')
    phi_pi      ${\phi_\pi}$            (long_name='inflation feedback')
    phi_y       ${\phi_y}$              (long_name='output feedback')
    phi_u       ${\phi_u}$              (long_name='unemployment feedback')
    @#endif
    ;


% Calibration page 1609
betta=0.99;
siggma=2;
epsilon=6;          % 20 percent markup
Psi=50;             % slope of PC
xi=0.4;             % matching function according to Blanchard/Diamond 
rho=0.08;           % empirical separation rate
varsigma=0.5;       % bargaining power
b_w_target=0.5;     % target for unemployment insurance
lambda_par=0.6;     % degree of wage rigidity
g_share=0.25;   
rho_Z=0.95;
rho_G=0.9;

@#if Ramsey==0
phi_R=0.9;          % empirical value
phi_pi=5;
phi_y=0;
phi_u=0;
@#endif

model;
[name='Definition marginal utility (3)']
lambda=c^-siggma;
[name='Euler Equation, eq. (4)']
1/R=betta*(lambda(+1)/lambda)/pi(+1);
[name='Definition labor market tightness, below (5)']
theta=v/u;
[name='Matching function, (5)']
m=m_param*(u^xi)*(v^(1-xi));
[name='Meeting rate between firms and workers, below (5)']
q=m/v;
[name='Aggregation of production function, eq. (6)']
y_gross=exp(z)*n;
[name='TFP process below eq. (6)']
z=rho_Z*z(-1)+epsilon_z;
[name='LOM employment, eq. (7)']
n=(1-rho)*(n(-1)+v(-1)*q(-1));
[name='Definition unemployment, eq. (8)']
u=1-n;
[name='FOC labor input, eq. (13)']
mu=mc*exp(z)-w+betta*lambda(+1)/lambda*(1-rho)*mu(+1);
[name='FOC vacancy postings, eq. (14)']
kappa/q=betta*(lambda(+1)/lambda)*(1-rho)*mu(+1);
[name='FOC price setting, eq. (15)']
1-Psi*(pi-1)*pi+betta*(lambda(+1)/lambda)*(Psi*(pi(+1)-1)*pi(+1)*y_gross(+1)/y_gross)=(1-mc)*epsilon;
[name='Nash wage schedule, Eq. (24)']
w=lambda_par*(varsigma*(mc*exp(z)+theta*kappa)+(1-varsigma)*b)+(1-lambda_par)*steady_state(w);
[name='Market clearing condition,']
y_net=c+g;
y_net=y_gross-kappa*v-y_gross*(Psi/2)*(pi-1)^2;
@#if Ramsey==0
    [name='Taylor rule, eq. (25)']
    log(R/steady_state(R))=phi_R*log(R(-1)/steady_state(R))+(1-phi_R)*(phi_pi*log(pi/1)+phi_y*log(y_net/y_net(-1))+phi_u*log(u/steady_state(u)));
@#endif

[name='Level of government spending']
g=G_SS*exp(z_G);
[name='LOM government spending']
z_G=rho_G*z_G(-1)+epsilon_G;

log_y_net=log(y_net);
log_v=log(v);
log_w=log(w);
log_u=log(u);
log_theta=log(theta);
log_pi=log(pi);
end;


initval;
R=1/betta;
end;

steady_state_model;
@#if Ramsey==0
    R=1/betta;
@#endif
pi=R*betta;             %implied by interest rate chosen by planner (but essentially fixed in price adjustment costs)
z=0; 
z_G=0;
mc=(epsilon-1+Psi*(pi-1)*pi*(1-betta))/epsilon; %collapses to (epsilon-1)/epsilon in a zero inflation steady state
q=0.7;                  %matching target
theta=0.6/q;            %matching target
m_param=q*theta^xi;     %implied by matching target
u=1/(1+q*theta*(1-rho)/rho);
n=1-u;
v=rho*n/((1-rho)*q);
m=m_param*(u^xi)*(v^(1-xi));
w=((1-(1-varsigma)*b_w_target)/(q*betta*(1-rho)*varsigma*theta)+1/(1-betta*(1-rho)))^(-1)*mc*(varsigma/(q*betta*(1-rho)*varsigma*theta)+1/(1-betta*(1-rho)));
mu=(mc-w)/(1-betta*(1-rho));
kappa=mu*(q*betta*(1-rho));
% kappa=(w-varsigma*mc-(1-varsigma)*b_w_target)/(varsigma*theta);
y_gross=n;
y_net=y_gross-kappa*v;
G_SS=g_share*y_net;
g=G_SS;
c=y_net-g;
lambda=c^-siggma;
b=b_w_target*w;
log_y_net=log(y_net);
log_v=log(v);
log_w=log(w);
log_u=log(u);
log_theta=log(theta);
log_pi=log(pi);
end;

@#if Ramsey
ramsey_model(planner_discount=0.99, instruments=(R));
@#endif

resid; 
steady; 
check;

write_latex_steady_state_model;
write_latex_static_model;
write_latex_dynamic_model(write_equation_tags);
write_latex_parameter_table;
collect_latex_files;

shocks; 
var epsilon_z; 
stderr 0.008;   %empirical estimate 
var epsilon_G; 
stderr 0.008;   %empirical estimate
end;

@#if Ramsey
    planner_objective(log(c));
@#endif

shocks; 
var epsilon_z; 
stderr 1;   %empirical estimate 
var epsilon_G; 
stderr 1;   %empirical estimate
end;

stoch_simul(irf=60, order=1, pruning, graph_format=pdf) log_theta log_v log_u log_w log_pi log_y_net;

fprintf('\nSteady state b/w: %4.3f\n',M_.params(strmatch('b',M_.param_names,'exact'))/oo_.dr.ys(strmatch('w',M_.endo_names,'exact')));
fprintf('Steady state theta*q: %4.3f\n',oo_.dr.ys(strmatch('theta',M_.endo_names,'exact'))*oo_.dr.ys(strmatch('q',M_.endo_names,'exact')));
fprintf('Steady state G/Y: %4.3f\n',oo_.dr.ys(strmatch('g',M_.endo_names,'exact'))/oo_.dr.ys(strmatch('y_net',M_.endo_names,'exact')));
fprintf('Steady state q: %4.3f\n',oo_.dr.ys(strmatch('q',M_.endo_names,'exact')));
fprintf('Steady state u: %4.3f\n',oo_.dr.ys(strmatch('u',M_.endo_names,'exact')));