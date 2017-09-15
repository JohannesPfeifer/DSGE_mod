/*
 * This file replicates some of the results in Jermann (1998): Asset pricing in production economies,
 * Journal of Monetary Economics, 41, pp. 257-275
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5
 *
 * Notes:
 *  - While Jermann uses log-linear asset pricing, this implementation relies on a 
 *      second order approximation
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model. 
 */

/*
 * Copyright (C) 2016-17 Johannes Pfeifer
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

var c           $C$             (long_name='consumption')
    k           $K$             (long_name='capital stock')
    invest      $I$             (long_name='investment')
    z           $z$             (long_name='productivity')
    lambda      $\lambda$       (long_name='marginal utility')
    q           $q$             (long_name='Tobins marginal q')
    w           $W$             (long_name='real wage')
    d           $d$             (long_name='dividends')
    V_k         ${V^k}$         (long_name='price of capital stock')
    V_b         ${V^b}$         (long_name='price of a consol')
    y           $Y$             (long_name='output')
    r_k         ${r^K}$         (long_name='return to capital')
    r_f         ${r^f}$         (long_name='risk-free interest rate')
    c_growth    ${\Delta C}$    (long_name='growth rate of consumption')
    y_growth    ${\Delta Y}$    (long_name='growth rate of output')
    i_growth    ${\Delta I}$    (long_name='growth rate of investment')
    log_d       ${\log(d)}$     (long_name='log dividend')
    log_V_k     ${\log(V^k)}$   (long_name='log price of capital stock')
    log_V_b     ${\log(V^b)}$   (long_name='log price of a consol')
    log_r_f     ${\log(r^f)}$   (long_name='log risk-free rate')
    log_r_k     ${\log(r^k)}$   (long_name='log return to capital')
    log_lambda  ${\log(\lambda)}$   (long_name='log marginal utility')
    SDF         ${\beta\frac{U_{C,t+1}}{U_{C,t}}}$  (long_name='stochastic discount factor')
    rf_ann      ${rf^{ann,net}}$    (long_name='annualized net risk free rate')
    rk_ann      ${rk^{ann,net}}$    (long_name='annualized net return to capital')
    rp_ann      ${ERP^{ann}}$   (long_name='annualized equity premium')
    equity_premium  ${ERP}$     (long_name='equity premium')
        ;

varexo e        ${\beta^*}$     (long_name='productivity shock')
        ;

predetermined_variables k;

parameters 
    betastar    ${\beta^*}$     (long_name='discount factor')
    delta       ${\delta}$      (long_name='depreciation rate')
    alpha       ${\alpha}$      (long_name='capital share')
    sigma       ${\sigma}$      (long_name='standard deviation productivity shock')
    rho         ${\rho}$        (long_name='autocorrelation productivity')
    gamma       ${\gamma}$      (long_name='long-run growth rate')
    tau         ${\tau}$        (long_name='risk aversion')
    h           ${h}$           (long_name='habit parameter (originally called alpha in the paper)')
    xi          ${\xi}$         (long_name='elasticity of investment-capital ratio wrt to q')
    const       ${c}$           (long_name='constant in in investment adjustment costs')
    a           ${a}$           (long_name='Parameter a in investment adjustment costs')
    b           ${b}$           (long_name='Parameter b in investment adjustment costs')
    i_k         ${\frac{I}{K}}$ (long_name='steady state investment to capital ratio')
    ;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
gamma   = 1.005;
alpha   = 0.36;
delta   = 0.025;
tau     = 5;
h       = 0.82;
betastar= gamma/1.011138;
sigma   = (0.0064/(1-alpha));
rho     = 0.99;
xi      = .23;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;
[name='1. Marginal utility']
lambda=(c-h/gamma*c(-1))^(-tau)-betastar*h/gamma*(c(+1)-h/gamma*c)^(-tau);
[name='2. Euler equation stocks']
lambda*V_k = betastar*lambda(+1)*(V_k(+1)+d(+1));
[name='Definition dividends']
d=exp(z)*k^alpha-w-invest;
[name='Real wage']
w=(1-alpha)*exp(z)*k^alpha;
[name='Resource constraint']
c+invest = exp(z)*k^alpha;
[name='LOM capital']
gamma*k(+1) = (1-delta)*k+(b/(1-a)*(invest/k)^(1-a)+const)*k;
[name='FOC capital']
lambda*q*gamma=betastar*lambda(+1)*(alpha*exp(z(+1))*k(+1)^(alpha-1)
                                    +q(+1)*(1-delta+const+b*a/(1-a)*(invest(+1)/k(+1))^(1-a)));
[name='FOC investment']
1=q*b*(invest/k)^-a;
[name='LOM technology']
z = rho*z(-1)+e;
[name='Production function']
y=exp(z)*k^alpha;
[name='Return to capital']
r_k=1/q*(alpha*exp(z(+1))*k(+1)^(alpha-1)+q(+1)*(1-delta+const+b*a/(1-a)*(invest(+1)/k(+1))^(1-a)));
[name='Definition stochastic discount factor']
SDF=betastar/gamma*lambda(+1)/lambda;
[name='Definition risk-free interest rate']
r_f= 1/SDF;
[name='log risk-free interest rate']
log_r_f=log(r_f); 
[name='log return to capital']
log_r_k=log(r_k); 

[name='log return to capital']
equity_premium = r_k - r_f;
[name='Annualized net risk-free rate']
rf_ann=4*(r_f-1); 
[name='Annualized net return to capital']
rk_ann=4*(r_k-1); 
[name='log return to capital']
rp_ann=4*equity_premium;

[name='Pricing equation for a consol']
lambda*V_b = betastar*lambda(+1)*(V_b(+1)+1/betastar-1);
[name='Definition growth rate of output']
y_growth=log(y)-log(y(-1));
[name='Definition growth rate of investment']
i_growth=log(invest)-log(invest(-1));
[name='Definition growth rate of consumption']
c_growth=log(c)-log(c(-1));

[name='Definition log price of capital']
log_V_k=log(V_k);
[name='Definition log price of consol']
log_V_b=log(V_b);
[name='Definition log dividend']
log_d=log(d);
[name='Definition log marginal utility']
log_lambda=log(lambda);
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

steady_state_model;
  a= 1/xi;//1/xi;  //used in functional form
  //beta=betastar*gamma^(1-tau);
  i_k=1-1/gamma*(1-delta);
  b=i_k^a;
  const=gamma*i_k-b/(1-a)*i_k^(1-a);
  k=((gamma/betastar-(1-delta+const + b*a/(1-a)*i_k^(1-a)))/alpha)^(1/(alpha-1));
  q=1;
  z = 0;
  invest=i_k*k;
  w=(1-alpha)*k^alpha;
  y=k^alpha;
  y_growth=0;
  d=k^alpha-w-invest;
  c=w+d;
  lambda=(c*(1-h/gamma))^(-tau)*(1-betastar*h/gamma);
  r_k=(alpha*k^(alpha-1)+q*(1-delta+const+b*a/(1-a)*(invest/k)^(1-a)));
  r_f= 1/(betastar/gamma);
  equity_premium = r_k - r_f;
  rf_ann=4*(r_f-1); 
  rk_ann=4*(r_k-1); 
  rp_ann=4*equity_premium;
  V_k=d/(1/betastar-1);
  V_b=1;
  log_V_k=log(V_k);
  log_V_b=log(V_b);
  log_d=log(d);
  log_r_f=log(r_f); 
  log_r_k=log(r_k); 
  log_lambda=log(lambda);
  m1   = (betastar/gamma);
  rf1  = 1/m1;
  SDF=betastar/gamma;
end;

write_latex_dynamic_model;

shocks;
var e = sigma^2; //picked to set output volatility to 1
end;

steady;
stoch_simul(order = 1,irf=0) y_growth c_growth i_growth ; 


stoch_simul(order = 2,irf=0,periods=50000) V_k d r_f r_k;

E_r_f=mean(exp(log_r_f)-1)*400

R=gamma*(exp(log_V_k(2:end))+exp(log_d(2:end)))./exp(log_V_k(1:end-1));
E_r_k=(mean(R)-1)*400

R_b=gamma*(exp(log_V_b(2:end))+1/betastar-1)./exp(log_V_b(1:end-1));
E_r_b=(mean(R_b)-1)*400

E_r_k-E_r_b