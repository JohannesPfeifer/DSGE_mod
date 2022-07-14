/*
 * This file replicates the IRFs in figure 3 of the RBC model with a 
 * constraint on investment illustrated in:
 * Guerrieri and Iacoviello (2015): "OccBin: A toolkit for solving dynamic 
 *   models with occastionally binding constraints easily", Journal of
 *   Monetary Economics 70, pp.22-38.
 *
 * Notes:
 *  - This mod file uses Dynare's occbin toolkit which is based on the one
 *    by Guerrieri and Iacoviello (2015)
 *  - The consumption Euler equation (11) in the paper has a typo at the 
 *    end as it should be \lambda_{t+1}; this is correct in the Appendix
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
  a        $A$            (long_name='technology')
  c        $C$            (long_name='consumption')
  iv       $I$            (long_name='investment')
  k        $K$            (long_name='capital')
  lam      $\lambda$      (long_name='Lagrange multiplier on investment constraint')  
  chat     ${\hat{c}}$    (long_name='consumption (percent dev from ss)')
  ivhat    ${\hat{i}}$    (long_name='investment (percent dev from ss)')
  khat     ${\hat{k}}$    (long_name='capital (percent dev from ss)')  
;
 
varexo
  epsi     $\epsi$        (long_name='negative TFP shock')
;

parameters
  ALPHA    $\alpha$       (long_name='bias towards capital in production')
  DELTA    $\delta$       (long_name='depreciation rate')
  BETA     $\beta$        (long_name='discount factor')
  GAMMA    $\gamma$       (long_name='risk aversion')
  RHO      $\rho$         (long_name='persistence technology')
  PHI      $\phi$         (long_name='fraction parameter for lower bound on investment')
;

ALPHA = 0.33;
DELTA = 0.1;
BETA  = 0.96;
GAMMA = 2;
RHO   = 0.9;
PHI   = 0.975;

model;
[name='consumption Euler equation, eq. (11), eq. (A.1)']
c^(-GAMMA) - lam = BETA*( c(+1)^(-GAMMA)*(1-DELTA+ALPHA*a(+1)*k^(ALPHA-1)) -(1-DELTA)*lam(+1) ); // note that there is a small typo in the paper, but not in the appendix
[name='resource constraint, eq. (7)']
c + iv = a*k(-1)^(ALPHA);
[name='capital accumulation, eq. (8)']
k = (1-DELTA)*k(-1) + iv;
[name='law of motion for technology, eq. (10)']
log(a) = RHO*log(a(-1)) + epsi;
[name='reporting: investment percent deviation from steady-state']
ivhat = 100*(iv/steady_state(iv)-1);
[name='reporting: consumption percent deviation from steady-state']
chat = 100*(c/steady_state(c)-1);
[name='reporting: capital percent deviation from steady-state']
khat = 100*(k/steady_state(k)-1);

[name='investment constraint',relax='irr']
lam = 0;
[name='investment constraint',bind='irr']
iv = PHI*steady_state(iv);
end;

steady_state_model;
a   = 1;
k   = ((1/BETA-1+DELTA)/ALPHA)^(1/(ALPHA-1));
c   = -DELTA*k + k^ALPHA;
iv  = DELTA*k;
lam = 0;
khat  = 0;
ivhat = 0;
chat  = 0;
end;
steady;

occbin_constraints;
name 'irr'; bind iv < PHI*steady_state(iv); relax lam <= 0;
end;

% replicate left column of Fig. 3
shocks(surprise,overwrite);
var epsi; periods 1; values -0.04;
end;
occbin_setup;
occbin_solver(simul_periods=50,simul_check_ahead_periods=100);
occbin_graph ivhat chat khat;

% replicate right column of Fig. 3
shocks(surprise,overwrite);
var epsi; periods 1; values 0.04;
end;
occbin_setup;
occbin_solver(simul_periods=100);
occbin_graph ivhat chat khat;