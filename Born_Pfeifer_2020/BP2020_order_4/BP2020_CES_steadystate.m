function [ys,params,check] =BP2020_CES_steadystate(ys,exe,M_,options_)
% function [ys,params,check] = BP2020_CES_steadystate(ys,exe,M_,options_)
% computes the steady state for the NK_baseline.mod and uses a numerical
% solver to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)

% Copyright (C) 2013-2020 Benjamin Born and Johannes Pfeifer
%
%  This is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  It is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  For a copy of the GNU General Public License,
%  see <http://www.gnu.org/licenses/>.

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
    paramname = deblank(M_.param_names{ii,:});
    eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;


%% Enter model equations here
theta_v =(1-siggma)/(1-1/chi); %definition utility function parameter
if separate_firm_risk_aversion_indicator
    theta_v_firm =(1-sigma_firm)/(1-1/chi_firm); %definition utility function parameter
end
V=1;
Y=1;
M=betta;
M_firm=M;
E_t_V_tp1_1_minus_sigma=1;
Xi = (theta_p-1)/theta_p;
mu=1/Xi;
Pi=Pi_bar;
R=Pi_bar/betta;
q=1;
u=1;
R_K=1/betta-(1-delta);
delta_1=R_K;
S=1;
Z=Z_bar;
sigma_z=sigma_z_bar;
sigma_G=sigma_G_bar;
G=1;

K=(1-labor_share)/R_K;
options=optimset('TolFun',10e-12,'TolX',10e-12,'Display','off');
%N=fsolve(@(N) (1-alppha)*(1-N-siggma*N)/(siggma*((2*N-1)/(1-N))*N)+delta*K-1,0.4)
if Iso_elastic_utility_indicator==0
    if ~steady_state_markup_indicator %0 markup due to infinite elasticity
        [N,fval,exitflag,output]=fsolve(@(N) (1-tau_l)/(1+tau_c)*labor_share/N*(1-N)*(theta_v/(1-siggma)*(1-Frisch_target*(1-(1-siggma)/theta_v)*N/(1-N)))/...
            (1-(theta_v/(1-siggma)*(1-Frisch_target*(1-(1-siggma)/theta_v)*N/(1-N))))...
            +delta*K-(1-gshare)*1,0.2,options);
    else %regular markup
        [N,fval,exitflag,output]=fsolve(@(N) ((theta_w-1)/theta_w)*(1-tau_l)/(1+tau_c)*labor_share/N*(1-N)*(theta_v/(1-siggma)*(1-Frisch_target*(1-(1-siggma)/theta_v)*N/(1-N)))/...
            (1-(theta_v/(1-siggma)*(1-Frisch_target*(1-(1-siggma)/theta_v)*N/(1-N))))...
            +delta*K-(1-gshare)*1,0.2,options);
    end
    if exitflag<1
        check=1;
        return
    end
    eta=(theta_v/(1-siggma)*(1-Frisch_target*(1-(1-siggma)/theta_v)*N/(1-N)));
    Frisch=(1-eta*(1-siggma)/theta_v)/(1-(1-siggma)/theta_v)*(1-N)/N;
    if abs(Frisch-Frisch_target)>10e-5
        error('Wrong calibration. Frisch elasticity is not at target')
    end
else
    N=1/3;
end
N_overhead=N_overhead_share*N;

if CES_indicator
    alppha=(1-labor_share)*(N-N_overhead)^psii_CES/(1-N_overhead_share)/(labor_share*K^psii_CES+(1-labor_share)*(N-N_overhead)^psii_CES/(1-N_overhead_share));
    PF_normalization=((alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES)* ...
        Xi*(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES*1/(1-N_overhead_share))/(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES))^(-1);
    Phi=PF_normalization*(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES)* ...
        (1-Xi*((alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES*1/(1-N_overhead_share))/(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)));
    MPN=PF_normalization*(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES-1)*(1-alppha)*(Z*(N-N_overhead))^psii_CES/(N-N_overhead);
    MPK=PF_normalization*(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES-1)*alppha*K^psii_CES/K;
    W=Xi*MPN;
    if abs((PF_normalization*(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES)-Phi)-1)>1e-10 ... output normalization to 1 not successful
            || abs(W*N/Y-labor_share)>1e-10 ...
            || abs(R_K*K/Y-(1-labor_share))>1e-10
        error('calibration wrong')
    end
else
    alppha=(1-labor_share)/(1-N_overhead_share)/(labor_share+(1-labor_share)/(1-N_overhead_share));
    PF_normalization=(1/(alppha*Xi*(N-N_overhead)^(1-alppha)*K^alppha+(1-alppha)*Xi*(N-N_overhead)^(1-alppha)*K^alppha*1/(1-N_overhead_share)));
    MPN=(1-alppha)*PF_normalization*K^alppha*(Z*(N-N_overhead))^(-alppha)*Z;
    MPK=alppha*PF_normalization*K^(alppha-1)*(Z*(N-N_overhead))^(1-alppha);
    Phi=(1-Xi*(alppha+(1-alppha)*1/(1-N_overhead_share)))*PF_normalization*(N-N_overhead)^(1-alppha)*K^alppha;
    Y=PF_normalization*K^alppha*(Z*(N-N_overhead))^(1-alppha)-Phi;
    W=Xi*MPN;
    if abs(Y-1)>1e-10 ... output normalization to 1 not successful
            || abs(W*N/Y-labor_share)>1e-10 ...
            || abs(R_K*K/Y-(1-labor_share))>1e-10
        error('calibration wrong')
    end
    
end

Y_net=Y; %net output including price and wage adjustment costs, but not fixed costs
D=Y-N*W-delta*K;
I=delta*K;
Gbar=gshare*Y;
C=Y-I-Gbar;

if Iso_elastic_utility_indicator==0
    V_normalization=(1-betta)/((C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v));

    U_C=V^(1-(1-siggma)/theta_v)*eta*V_normalization*(C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v)/C;
    U_N=-(V_normalization*(C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v)+betta*(E_t_V_tp1_1_minus_sigma)^(1/theta_v))^(theta_v/(1-siggma)-1)*V_normalization*(C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v)*(1-eta)/(1-N);
else
    eta=((1-h)*C)^-siggma/(N^Frisch_target)*W*(1-tau_l)/(1+tau_c)*(theta_w-1)/theta_w;

    V_normalization=(1-betta)/((C*(1-h))^(1-siggma)/(1-siggma)-eta*N^(1+Frisch_target)/(1+Frisch_target));

    U_C=V_normalization*(C*(1-h))^(-siggma);
    U_N=-V_normalization*eta*N^(Frisch_target);   
end

MRS=-U_N/U_C; %equal to (1-eta)/eta*C/(1-N);

lambda=U_C/(1+tau_c);

T=tau_c*C+tau_l*W*N-Gbar;

firm_wedge=log(MPN/W);
if steady_state_markup_indicator && ~overhead_labor_indicator %0 markup due to infinite elasticity
    if abs(firm_wedge-log(theta_p/(theta_p-1))) >1e-10;
        error('Firm Wedge is wrong')
    end
end
household_wedge=-log(MRS*((1+tau_c)/(1-tau_l))/W);
if steady_state_markup_indicator %0 markup due to infinite elasticity
    if abs(household_wedge-log(theta_w/(theta_w-1)))>1e-10
        error('Household Wedge is wrong')
    end
end

Y_hp_cyc=0;
log_Y=log(Y);
log_C=log(C);
log_I=log(I);
log_mu=log(mu);
log_N=log(N);
log_W =log(W);
log_K =log(K);
pi_annualized=1+(4*(Pi-1));
R_annualized=1+(4*(R-1));
R_K_annualized=4*R_K;
log_pi_annualized=log(pi_annualized);
log_R_annualized=log(R_annualized);
log_R_K_annualized=log(R_K_annualized);

%% end own model equations
for iter = 1:length(M_.params) %update parameters set in the file
    eval([ 'params(' num2str(iter) ',1) = ' M_.param_names{iter,:} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
    varname = deblank(M_.endo_names{ii,:});
    eval(['ys(' int2str(ii) ') = ' varname ';']);
end
