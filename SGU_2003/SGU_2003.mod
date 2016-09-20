/*
 * This file replicates the model studied in:
 * Schmitt-Grohé, Stephanie and Uribe, Martín (2003): "Closing small open economy 
 * models", Journal of International Economics, 61, pp. 163-185.
 * 
 * It provides a replication code for the main results of the original paper  
 * for Models 1 (debt elastic interest rate) to 4 (portfolio holding costs). 
 *
 * Notes:
 * - This file has been tested with Dynare 4.4.3
 * - SGU on page 169 define the trade balance as 1-(c+i)/y. This in principle neglects the resource costs 
 *      imposed by capital adjustment costs (thanks to Sidney Caetano for pointing this out). The present mod-
 *      file therefore uses the correct definition     
 *          tb_y = 1-(c+i+(phi/2)*(k-k(-1))^2)/y;
 *      Note that this change does not affect any of the results presented in the paper, because they were obtained at first order
 *      and the first order derivative of this term is 0.
 *   
 * This implementation was written by Johannes Pfeifer. If you spot any mistakes, 
 * email me at jpfeifer@gmx.de.
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2015 Johannes Pfeifer
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

% Endogenous discount factor
@#define model1 =1
% Endogenous discount factor without internalization
@#define model1a =0
%debt elastic interest rate premium
@#define model2 =0
%portfolio holding costs
@#define model3 =0 
%Complete markets
@#define model4 =0 

var  c h y i k a lambda ${\lambda}$ util;  
varexo e;                                    
                                             
parameters  gamma ${\gamma}$
            omega ${\omega}$
            rho ${\gamma}$
            sigma_tfp ${\sigma_{a}}$
            delta ${\delta}$
            psi_1 ${\psi_1}$
            psi_2 ${\psi_2}$
            alpha ${\alpha}$
            phi ${\phi}$
            psi_3 ${\psi_3}$
            psi_4 ${\psi_4}$
            r_bar ${\bar r}$
            d_bar ${\bar d}$;
            
%Table 1
gamma  = 2; %risk aversion
omega  = 1.455; %Frisch-elasticity parameter
psi_1  = 0; %set in steady state %elasticity discount factor w.r.t. to arguments of utility function
alpha  = 0.32; %labor share
phi    = 0.028; %capital adjustment cost parameter
r_bar    = 0.04; %world interest rate		
delta  = 0.1; %depreciation rate
rho    = 0.42; %autocorrelation TFP 
sigma_tfp = 0.0129; %standard deviation TFP

%Table 2
psi_2    = 0.000742;
d_bar  = 0.7442;


psi_3  = 0.00074; 
psi_4  = 0; %set in steady state; parameter complete markets case


@#if model1 == 1
var d tb_y, ca_y, r beta_fun, eta;   

model;
    //1. Eq. (4), Evolution of debt
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2;
    //2. Eq. (5), Production function
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    //3. Eq. (6), Law of motion for capital
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 

    //4. Eq. (8), Euler equation
    exp(lambda)= beta_fun*(1+exp(r))*exp(lambda(+1)); 
    //5. Eq. (9), Definition marginal utility
    exp(lambda)=(exp(c)-((exp(h)^omega)/omega))^(-gamma)-exp(eta)*(-psi_1*(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1-1));  
    //6. Eq. (10), Law of motion Lagrange mulitplier on discount factor equation
    exp(eta)=-util(+1)+exp(eta(+1))*beta_fun(+1);
    //7. Eq. (11), Labor FOC
    ((exp(c)-(exp(h)^omega)/omega)^(-gamma))*(exp(h)^(omega-1)) + 
        exp(eta)*(-psi_1*(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1-1)*(-exp(h)^(omega-1))) = exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    //8. Eq. (12), Investment FOC
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta_fun*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    //9. Eq. (14), Law of motion for TFP
    a = rho*a(-1)+sigma_tfp*e; 

    //10. Definition endogenous discount factor, p. 168
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    //9. Eq. (23), country interest rate 
    exp(r) = r_bar;

    //11. p. 169, Definition of trade balance to ouput ratio
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    ca_y = (1/exp(y))*(d(-1)-d);                                   
end;

steady_state_model;
    r     = log(r_bar);
    d     = d_bar;
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = log(exp(y)-exp(i)-r_bar*d);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    psi_1=-log(1/(1+r_bar))/(log((1+exp(c)-omega^(-1)*exp(h)^omega)));
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
    eta=log(-util/(1-beta_fun));
    lambda=log((exp(c)-((exp(h)^omega)/omega))^(-gamma)-exp(eta)*(-psi_1*(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1-1)));
    a     = 0;
    ca_y    = 0;
end;

@# endif

@#if model1a == 1
var d tb_y, ca_y, r beta_fun;   

model;
    //1. Eq. (4), Evolution of debt
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2;
    //2. Eq. (5), Production function
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    //3. Eq. (6), Law of motion for capital
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 

    //4. Eq. (17), Euler equation
    exp(lambda)= beta_fun*(1+exp(r))*exp(lambda(+1)); 
    //5. Eq. (18), Definition marginal utility
    exp(lambda)=(exp(c)-((exp(h)^omega)/omega))^(-gamma);  
    //6. Eq. (11), Labor FOC
    ((exp(c)-(exp(h)^omega)/omega)^(-gamma))*(exp(h)^(omega-1))= exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    //7. Eq. (12), Investment FOC
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta_fun*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    //8. Eq. (14), Law of motion for TFP
    a = rho*a(-1)+sigma_tfp*e; 

    //9. Definition endogenous discount factor, p. 168
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    //10. Eq. (23), country interest rate 
    exp(r) = r_bar;

    //12. p. 169, Definition of trade balance to ouput ratio
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    //13. Definition of current account to ouput ratio
    ca_y = (1/exp(y))*(d(-1)-d);                                   
end;

steady_state_model;
    r     = log(r_bar);        
    d     = d_bar;
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = log(exp(y)-exp(i)-r_bar*d);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    psi_1=-log(1/(1+r_bar))/(log((1+exp(c)-omega^(-1)*exp(h)^omega)));
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
    a     = 0;
    ca_y    = 0;
end;

@# endif

@#if model2 == 1
var d tb_y, ca_y, r riskpremium;
parameters beta ${\beta}$;

model;
    //1. Eq. (4), Evolution of debt
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2;
    //2. Eq. (5), Production function
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    //3. Eq. (6), Law of motion for capital
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 

    //4. Eq. (24), Euler equation
    exp(lambda)= beta*(1+exp(r))*exp(lambda(+1)); 
    //5. Eq. (25), Definition marginal utility
    (exp(c)-((exp(h)^omega)/omega))^(-gamma)   = exp(lambda);  
    //6. Eq. (26), Labor FOC
    ((exp(c)-((exp(h)^omega)/omega))^(-gamma))*(exp(h)^(omega-1))  = exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    //7. Eq. (27), Investment FOC
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    //8. Eq. (14), Law of motion for TFP
    a = rho*a(-1)+sigma_tfp*e; 
    //9. Eq. (23), country interest rate 
    exp(r) = r_bar+riskpremium;
    //10. p. 171 below Eq. (28), definition risk premia
    riskpremium = psi_2*(exp(d-d_bar)-1);

    //11. p. 169, Definition of trade balance to ouput ratio
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    ca_y = (1/exp(y))*(d(-1)-d);                                   
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
end;

steady_state_model;
    beta  = 1/(1+r_bar);
    r     = log((1-beta)/beta);
    d     = d_bar;
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = log(exp(y)-exp(i)-r_bar*d);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
    a     = 0;
    ca_y    = 0;
    riskpremium = 0;
end;

@# endif



@#if model3 == 1
var d tb_y, ca_y, r;
parameters beta ${\beta}$;

model;
    //1. Eq. (29), Evolution of debt
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2+psi_3*(d-d_bar)^2;
    //2. Eq. (5), Production function
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    //3. Eq. (6), Law of motion for capital
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 

    //1. Eq. (30),  Euler equation
    exp(lambda)*(1-psi_3*(d-d_bar))= beta*(1+exp(r))*exp(lambda(+1)); 
    //5. Eq. (25), Definition marginal utility
    (exp(c)-((exp(h)^omega)/omega))^(-gamma)   = exp(lambda);  
    //6. Eq. (26), Labor FOC
    ((exp(c)-((exp(h)^omega)/omega))^(-gamma))*(exp(h)^(omega-1))  = exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    //7. Eq. (27), Investment FOC
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    //8. Eq. (14), Law of motion for TFP
    a = rho*a(-1)+sigma_tfp*e; 
    //9. Eq. (13), country interest rate 
    exp(r) = r_bar;

    //11. p. 169, Definition of trade balance to ouput ratio
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    ca_y = (1/exp(y))*(d(-1)-d); 
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
end;

steady_state_model;
    beta  = 1/(1+r_bar);
    r     = log((1-beta)/beta);
    d     = d_bar;
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = log(exp(y)-exp(i)-r_bar*d);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
    a     = 0;
    ca_y    = 0;
end;

@# endif

@#if model4 == 1
    var tb_y;
parameters beta ${\beta}$;

model;
    //1. Eq. (5), Production function
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    //2. Eq. (6), Law of motion for capital
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 
    //3. Eq. (25), Definition marginal utility
    (exp(c)-((exp(h)^omega)/omega))^(-gamma)   = exp(lambda);  
    //4. Eq. (26), Labor FOC
    ((exp(c)-((exp(h)^omega)/omega))^(-gamma))*(exp(h)^(omega-1))  = exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    //5. Eq. (27), Investment FOC
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    //6. Eq. (35),  Euler equation
    exp(lambda)= psi_4; 
    //7. Eq. (14), Law of motion for TFP
    a = rho*a(-1)+sigma_tfp*e; 
    //8. p. 169, Definition of trade balance to ouput ratio
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
end;

steady_state_model;
    beta  = 1/(1+r_bar);
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = 0.110602; %from incomplete markets case
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
    psi_4=exp(lambda);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    a     = 0;
end;

@# endif


resid(1);

check;
steady; 


shocks;
    var e; stderr 1;
end;

stoch_simul(order=1, irf=0);

//Report results from Table 3
y_pos=strmatch('y',M_.endo_names,'exact');
c_pos=strmatch('c',M_.endo_names,'exact');
i_pos=strmatch('i',M_.endo_names,'exact');
h_pos=strmatch('h',M_.endo_names,'exact');
tb_y_pos=strmatch('tb_y',M_.endo_names,'exact');
ca_y_pos=strmatch('ca_y',M_.endo_names,'exact');


fprintf('\nstd(y):              \t %2.1f \n',sqrt(oo_.var(y_pos,y_pos))*100)
fprintf('std(c):                \t %2.1f \n',sqrt(oo_.var(c_pos,c_pos))*100)
fprintf('std(i):                \t %2.1f \n',sqrt(oo_.var(i_pos,i_pos))*100)
fprintf('std(h):                \t %2.1f \n',sqrt(oo_.var(h_pos,h_pos))*100)
fprintf('std(tb/y):             \t %2.1f \n',sqrt(oo_.var(tb_y_pos,tb_y_pos))*100)
if ~isempty(ca_y_pos)
fprintf('std(ca/y):             \t %2.1f \n',sqrt(oo_.var(ca_y_pos,ca_y_pos))*100)
else %complete markets case
fprintf('std(ca/y):             \t %2.1f \n',sqrt(oo_.var(ca_y_pos,ca_y_pos))*100)
end
fprintf('corr(y_t,y_t-1):       \t %3.2f \n',oo_.autocorr{1}(y_pos,y_pos))
fprintf('corr(c_t,c_t-1):       \t %3.2f \n',oo_.autocorr{1}(c_pos,c_pos))
fprintf('corr(i_t,i_t-1):       \t %4.3f \n',oo_.autocorr{1}(i_pos,i_pos))
fprintf('corr(h_t,h_t-1):       \t %3.2f \n',oo_.autocorr{1}(h_pos,h_pos))
fprintf('corr(tb/y_t,tb/y_t-1): \t %3.2f \n',oo_.autocorr{1}(tb_y_pos,tb_y_pos))
if ~isempty(ca_y_pos)
fprintf('corr(ca/y_t,ca/y_t-1): \t %3.2f \n',oo_.autocorr{1}(ca_y_pos,ca_y_pos))
else %complete markets case
fprintf('corr(ca/y_t,ca/y_t-1): \t %3.2f \n',NaN)
end
fprintf('corr(c_t,y_t):         \t %3.2f \n',oo_.gamma_y{1}(c_pos,y_pos)/sqrt(oo_.var(c_pos,c_pos)*oo_.var(y_pos,y_pos)))
fprintf('corr(i_t,y_t):         \t %3.2f \n',oo_.gamma_y{1}(i_pos,y_pos)/sqrt(oo_.var(i_pos,i_pos)*oo_.var(y_pos,y_pos)))
fprintf('corr(h_t,y_t):         \t %2.1f \n',oo_.gamma_y{1}(h_pos,y_pos)/sqrt(oo_.var(h_pos,h_pos)*oo_.var(y_pos,y_pos)))
fprintf('corr(tb/y_t,y_t):      \t %4.3f \n',oo_.gamma_y{1}(tb_y_pos,y_pos)/sqrt(oo_.var(tb_y_pos,tb_y_pos)*oo_.var(y_pos,y_pos)))
if ~isempty(ca_y_pos)
fprintf('corr(ca/y_t,y_t):      \t %4.3f \n',oo_.gamma_y{1}(ca_y_pos,y_pos)/sqrt(oo_.var(ca_y_pos,ca_y_pos)*oo_.var(y_pos,y_pos)))
else %complete markets case
fprintf('corr(ca/y_t,y_t):      \t %4.3f \n',NaN)
end


// Generate IRFs in Figure 1
shocks;
    var e; stderr 1/sigma_tfp; %Normalize to unit shock
end;

stoch_simul(order=1, irf=10);
