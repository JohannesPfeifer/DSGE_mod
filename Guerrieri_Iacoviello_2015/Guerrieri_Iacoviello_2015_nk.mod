/*
 * This file replicates the IRFs in figure 5 of the NK model with a 
 * constraint on nominal interest rates illustrated in:
 * Guerrieri and Iacoviello (2015): "OccBin: A toolkit for solving dynamic 
 *   models with occastionally binding constraints easily", Journal of
 *   Monetary Economics 70, pp.22-38.
 *
 * Notes:
 *  - This mod file uses Dynare's occbin toolkit which is based on the one
 *    by Guerrieri and Iacoviello (2015)
 *  - The model equations are taken from the appendix as the equations in 
 *    the printed version have some typos, e.g. missing paranthesis in eq. (24).
 *    Moreover, in the original replication files the model equations were
 *    simplified to reduce the model size for the global solution method;
 *    however, for this replication this is not important so we stick with
 *    the full version as outlined in the appendix.
 *  - The authors do not provide a calibrated value for PSI. Looking at 
 *    the original replication files it looks like PSI is an endogenous
 *    parameter that makes sure that in steady-state y=1. This is implemented
 *    in the steady-state model block in this mod file.
 *  - The value for the shock is not clear.
 *    In the paper it is claimed that the size is equal to 4 standard 
 *    deviations and the standard deviation is given in table 3 as 0.005.
 *    However, a positive shock of size 4*0.005=0.02 brings \beta_6 up to 
 *    1.000514 (and not to 1.019 as claimed in Fig. 5) and a negative shock
 *    of this size brings \beta_6 down to 0.987486 (and not to 0.969 as 
 *    claimed in Fig. 5). In the replication files a value of 0.024 is used;
 *    however, the positive shock of this size brings \beta_6 up to 1.001817
 *    and a negative shock down to 0.986183. The printed figure, however,
 *    looks closest to a value of 5*0.005=0.025, which is chosen below.
 *    Of course one can change the size below.
 *  - The reporting variables are different in the replication files than
 *    in the published paper. This implementation focuses on the figure
 *    in the published paper.
 * 
 * This implementation was written by Willi Mutschler (@wmutschl).
 *
 * If you spot mistakes, email me at willi@mutschler.eu
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright © 2022 Willi Mutschler
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
 * see <https://www.gnu.org/licenses/>.
 */

var
  bet        ${\beta}$              (long_name='time-varying discount factor')
  c          ${C}$                  (long_name='consumption')
  y          ${Y}$                  (long_name='output')
  l          ${L}$                  (long_name='labor')
  w          ${w}$                  (long_name='wage')
  mc         ${mc}$                 (long_name='marginal cost')  
  r          ${R}$                  (long_name='nominal interest rate')
  g          ${G}$                  (long_name='government spending')
  pie        ${\Pi}$                (long_name='inflation')
  pie_star   ${\Pi^*}$              (long_name='reset price')
  x1         ${x1}$                 (long_name='aux. sum 1 recursive price setting')
  x2         ${x2}$                 (long_name='aux. sum 1 recursive price setting')
  v          ${v}$                  (long_name='price dispersion')
  pie_an     ${\pi^{ann}}$          (long_name='annualized inflation (in percentage point deviation from steady-state)')
  r_an       ${R^{ann}}$            (long_name='annualized nominal interest rate (in level percent deviation from ZLB)')
  yhat       ${\hat{y}}$            (long_name='output (in percentage deviation from steady-state)')
;

varexo
  epsi       ${\epsilon}$           (long_name='discount factor shock')
;

parameters
  BETA       ${\beta}$              (long_name='steady-state discount factor')
  THETA      ${\theta}$             (long_name='Calvo parameter')
  PHI_Y      ${\phi_y}$             (long_name='response to output, Mon. Pol. Rule')
  PHI_PI     ${\phi_\pi}$           (long_name='response to inflation, Mon. Pol. Rule')
  RHO        ${\rho}$               (long_name='persistence of discount rate shock')
  EPSILON    ${\varepsilon}$        (long_name='elasticity of substitution across goods')
  G_Y        ${s_g}$                (long_name='steady-state ratio of G/Y')
  PI         ${\Pi}$                (long_name='steady-state inflation')
  PHI        ${\phi}$               (long_name='labor supply elasticity')
  PSI        ${\psi}$               (long_name='utility weight leisure (endogenous)')
  ZLB        ${ZLB}$                (long_name='value of zero lower bound')
;

BETA    = 0.994;
THETA   = 0.90;
PHI_Y   = 0.25;
PHI_PI  = 2.5;
RHO     = 0.8;
EPSILON = 6;
G_Y     = 0.20;
PI      = 1.005;
PHI     = 1;
ZLB     = 1;
% note that PSI is calibrated in the steady_state model block such that y=1;

model;
[name='(A.2): consumption Euler equation']
c^(-1) = bet*c(+1)^(-1)*r/pie(+1);
[name='(A.3): marginal costs / labor demand']
mc = w;
[name='(A.4): labor supply']
w = PSI*l^PHI*c;
[name='(A.5): optimal price setting']
EPSILON*x1 = (EPSILON-1)*x2;
[name='(A.6): optimal price setting auxiliary recursion 1']
x1 = mc*y/c + THETA*bet*pie(+1)^EPSILON*x1(+1);
[name='(A.7): optimal price setting auxiliary recursion 2']
x2/pie_star = y/c + THETA*bet*pie(+1)^(EPSILON-1)*x2(+1)/pie_star(+1);
[name='(A.8): monetary policy rule',relax='zlb']
r = steady_state(r)*(pie/PI)^PHI_PI*(y/steady_state(y))^PHI_Y;
[name='(A.8): monetary policy rule',bind='zlb']
r = ZLB;
[name='(A.10): government spending']
g = G_Y*y;
[name='(A.11): law of motion for optimal reset price']
1 = THETA*pie^(EPSILON-1) + (1-THETA)*pie_star^(1-EPSILON);
[name='(A.12): law of motion for price dispersion']
v = THETA*pie^EPSILON*v(-1) + (1-THETA)*pie_star^(-EPSILON);
[name='(A.13): aggregate demand']
y = c + g;
[name='(A.14): aggregate supply']
y = l/v;
[name='(A.15): law of motion discount factor']
log(bet) = (1-RHO)*log(BETA) + RHO*log(bet(-1)) + epsi;

[name='reporting: annualized nominal interest rate (in level percent deviation from ZLB)']
r_an = 400*(r-ZLB);
[name='reporting: reporting: annualized inflation (in percentage point deviation from steady-state)']
pie_an = 400*(pie-steady_state(pie));
[name='reporting: output (in percentage deviation from steady-state)']
yhat = 100*(y/steady_state(y)-1);
end;


steady_state_model;
bet = BETA;
pie = PI;
pie_star = ((1-THETA*pie^(EPSILON-1))/(1-THETA))^(1/(1-EPSILON));
mc = (EPSILON-1)/EPSILON * pie_star * (1-bet*THETA*pie^EPSILON) / (1-bet*THETA*pie^(EPSILON-1)) ;
v = (1-THETA)*pie_star^(-EPSILON)/(1-THETA*pie^EPSILON);
r = pie/BETA;
w = mc;
y_c = 1/(1-G_Y);
x2 = pie_star*y_c / (1-THETA*bet*pie^(EPSILON-1));
x1 = (EPSILON-1)/EPSILON*x2;
% calibrate PSI such that y=1 in steady-state
y = 1;
l = v*y;
c = (1-G_Y)*y;
g = G_Y*y;
PSI = w/(l^PHI*c); % note that this is an endogenous parameter updated during steady-state computations
r_an = 400*(r-ZLB);
pie_an = 0;
yhat = 0;
end;
steady;

occbin_constraints;
name 'zlb'; bind r <=  ZLB; relax r > ZLB;
end;

% shock size: see notes on top
@#define SHOCKSIZE = 5*0.005

% replicate left column of Fig. 5 (increase in discount factor)
shocks(surprise,overwrite);
var epsi; periods 6; values @{SHOCKSIZE};
end;

occbin_setup(simul_periods=30,simul_check_ahead_periods=200);
occbin_solver;
% the figure only plots simulations for periods 6 to 15
oo_.occbin.simul.piecewise = oo_.occbin.simul.piecewise(6:15,:);
oo_.occbin.simul.linear = oo_.occbin.simul.linear(6:15,:);
occbin_graph r_an v pie_an yhat;
fprintf('positive shock of size %f brings beta up to %f\n',options_.occbin.simul.SHOCKS(6,1),oo_.occbin.simul.piecewise(6,ismember(M_.endo_names,'bet')));

% replicate left column of Fig. 5 (decrease in discount factor)
shocks(surprise,overwrite);
var epsi; periods 6; values -@{SHOCKSIZE};
end;
occbin_setup(simul_periods=30,simul_check_ahead_periods=200);
occbin_solver;
% the figure only plots simulations for periods 6 to 15
oo_.occbin.simul.piecewise = oo_.occbin.simul.piecewise(6:15,:);
oo_.occbin.simul.linear = oo_.occbin.simul.linear(6:15,:); % 
occbin_graph r_an v pie_an yhat;
fprintf('negative shock of size %f brings beta down to %f\n',options_.occbin.simul.SHOCKS(6,1),oo_.occbin.simul.piecewise(6,ismember(M_.endo_names,'bet')));

