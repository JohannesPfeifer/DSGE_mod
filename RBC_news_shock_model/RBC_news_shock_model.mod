/*
 * This file implements a simple RBC model with additively separable utility and TFP news calibrated to US data. It shows how 
 * to generate IRFs to a "pure" news shock where an 8 period anticipated news shock does not materialize at time 0. 
 * This is the type of policy experiment that is for example performed in Beaudry Portier (2004): An exploration  
 * into Pigou’s theory of cycles, Journal of Monetary Economics 51, pp. 1183–1216.
 *
 * This is done using Dynare's simult_ function, which can be used to simulate time series given the decision rules. This capability
 * can be useful if user's want to use a particular shock series (e.g. truncated normal, historical smoothed shocks etc.)
 * 
 * Note that having two exactly offsetting shocks from continuous distributions to keep the exogenous variable constant 
 * is a 0 probability event for the agents of the model. 
 * 
 * Moreover, this mod-file shows how to use Dynare's capacities to generate TeX-files of the model equations. If you want to see the model
 * equations belonging to this mod-file, run it using Dynare and then use a TeX-editor to compile the TeX-files generated.
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * This mod-file uses variable substitution to perform a log-linearization of the model. All variables,
 * except for the interest rate (which is already in percent), are put in exp() for this purpose.
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.

 */

/*
 * Copyright (C) 2013-15 Johannes Pfeifer
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


var y c k l z r w 
    invest ${i}$; //define that the variable invest should be called i in the TeX-file

varexo eps_z_news ${\varepsilon_z^{news}}$ //use out brackets to prevent problems with multiple sub- and superscripts
    eps_z_surprise ${\varepsilon_z^{surprise}}$;


parameters beta $\beta$
    psi $\psi$
    sigma $\sigma$
    delta $\delta$
    alpha $\alpha$
    rhoz $\rho_z$
    gammax $\gamma_x$
    n
    x 
    i_y 
    k_y;

% set parameter values
sigma=1;                % risk aversion
alpha= 0.33;            % capital share
i_y=0.25;               % investment-output ration
k_y=10.4;               % capital-output ratio
x=0.0055;               % technology growth (per capita output growth)
n=0.0027;               % population growth
rhoz=0.97;              %technology autocorrelation base on linearly detrended Solow residual

model;
//1. Euler equation
exp(c)^(-sigma)=beta/gammax*exp(c(+1))^(-sigma)* (alpha*exp(z(+1))*(exp(k)/exp(l(+1)))^(alpha-1)+(1-delta));
//2. Labor FOC
psi*exp(c)^sigma*1/(1-exp(l))=exp(w);
//3. Law of motion capital 
gammax*exp(k)=(1-delta)*exp(k(-1))+exp(invest);
//4. resource constraint
exp(y)=exp(invest)+exp(c);
//5. production function
exp(y)=exp(z)*exp(k(-1))^alpha*exp(l)^(1-alpha);
//6. real wage/firm FOC labor
exp(w)=(1-alpha)*exp(y)/exp(l);
//7. annualized real interst rate/firm FOC capital
r=4*alpha*exp(y)/exp(k(-1));
//8. exogenous TFP process
z=rhoz*z(-1)+eps_z_surprise + eps_z_news(-8);
end;

steady_state_model; %use steady_state_model-block to calibrate model and set parameters
    gammax=(1+n)*(1+x); 
    % set depreciation rate and discount factor to be compatible with steady state capital-output and investment-output ratio
    delta=i_y/k_y-x-n-n*x;  %deprecation rate
    beta=gammax/(alpha/k_y+(1-delta)); %discount factor
    % calibrate the model to steady state labor of 0.33, i.e. compute the corresponding steady state values
    % and the labor disutility parameter by hand; the steady state values are used later in the initval-block
    l_ss=0.33;
    k_ss=((1/beta*gammax-(1-delta))/alpha)^(1/(alpha-1))*l_ss;
    i_ss=(x+n+delta+n*x)*k_ss;
    y_ss=k_ss^alpha*l_ss^(1-alpha);
    c_ss=y_ss-i_ss;
    psi=(1-alpha)*(k_ss/l_ss)^alpha*(1-l_ss)/c_ss^sigma; %labor disutility parameter that sets labor to 0.33 in steady state
    invest=log(i_ss);
    w=log((1-alpha)*y_ss/l_ss);
    r=4*alpha*y_ss/k_ss;
    y=log(y_ss);
    k=log(k_ss);
    c=log(c_ss);
    l=log(l_ss);
    z=0;
end;

shocks;
    var eps_z_news=1; //8 period anticipated TFP news shock
    var eps_z_surprise=1; //TFP surprise shock
end;

//use Dynare capabilities to generate TeX-files of the dynamic and static model
write_latex_static_model;
write_latex_dynamic_model;

steady;
check;

stoch_simul(order=1,irf=40);


//************ The following line generate the IRFs to a "pure" TFP news shock, i.e. the anticipated shock is 
//************ counteracted by an opposite surprise shock when it is supposed to realize

//initialize IRF generation
initial_condition_states = repmat(oo_.dr.ys,1,M_.maximum_lag);
shock_matrix = zeros(options_.irf,M_.exo_nbr); %create shock matrix with number of time periods in columns

// set shocks for pure news 
shock_matrix(1,strmatch('eps_z_news',M_.exo_names,'exact')) = 1; %set news shock to 1 (use any shock size you want)
shock_matrix(1+8,strmatch('eps_z_surprise',M_.exo_names,'exact')) = -1; %8 periods later use counteracting shock of -1

y2 = simult_(initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,options_.irf); %deviation from steady state


// manually select variables for figure
figure
subplot(2,1,1)
plot(y_IRF(strmatch('y',M_.endo_names,'exact'),:)); % use strmatch to select values
title('Output');
subplot(2,1,2)
plot(y_IRF(strmatch('z',M_.endo_names,'exact'),:));
title('TFP');

// Automatically loop over variables for figure (may require different setting for subplots in larger models)
figure
for ii=1:M_.orig_endo_nbr
    subplot(3,3,ii)
    if max(abs(y_IRF(ii,:)))>1e-12 %get rid of numerical inaccuracies
        plot(y_IRF(ii,:));
    else
        plot(zeros(options_.irf,1));
    end
    title(deblank(M_.endo_names(ii,:)));
end

