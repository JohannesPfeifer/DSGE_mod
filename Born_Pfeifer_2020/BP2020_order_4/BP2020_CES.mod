/* 
 * Born/Pfeifer (2020) model with sticky wages and wages as well as overhead labor and CES production function
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.6.2 OR HIGHER
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model. 
*/

/*
 * Copyright (C) 2020 Johannes Pfeifer and Benjamin Born
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

@#ifndef Large_shocks 
    @#define Large_shocks = 0
@#endif

@#ifndef Sequential_shocks 
    @#define Sequential_shocks = 0
@#endif

@#ifndef Leduc_Liu_calibration 
    @#define Leduc_Liu_calibration = 0
@#endif

@#ifndef FV_et_al_calibration 
    @#define FV_et_al_calibration = 0
@#endif

@#ifndef sticky_prices 
    @#define sticky_prices = 0
@#endif

@#ifndef FVEtal_Rotemberg 
    @#define FVEtal_Rotemberg = 1
@#endif

@#ifndef sticky_wages
    @#define sticky_wages = 0
@#endif

@#ifndef HP_Taylor 
    @#define HP_Taylor = 1
@#endif

@#ifndef CES 
    @#define CES = 1
@#endif

@#ifndef overhead_labor 
    @#define overhead_labor = 1
@#endif
    
@#ifndef Adjustments_Costs 
    @#define Adjustments_Costs = 1
@#endif
    
@#ifndef CEE 
        @#define CEE = 1    
@#endif

@#ifndef capital_utilization 
    @#define capital_utilization = 0
@#endif

@#ifndef CRRA
    @#define CRRA = 1
@#endif

@#ifndef Iso_elastic_utility 
    @#define Iso_elastic_utility = 0
@#endif

@#ifndef Separate_Firm_Risk_Aversion
    @#define Separate_Firm_Risk_Aversion = 0
@#endif

var C       ${C}$                               (long_name='consumption')
    N       ${N}$                               (long_name='labor')
    W       ${\left({\frac{W}{P}}\right)}$      (long_name='real wage')
    V       ${V}$                               (long_name='Value function')
    E_t_V_tp1_1_minus_sigma ${\left({E_t V_{t+1}^{1-\sigma}}\right)}$ (long_name='Auxilary variable EZ')
    Y       ${Y}$                               (long_name='Output')
    K       ${K}$                               (long_name='capital')
    Z       ${Z}$                               (long_name='TFP')
    I       ${I}$                               (long_name='Investment')
    D       ${\frac{D}{P}}$                     (long_name='Cash flows of the firm')
    Pi      ${\Pi}$                             (long_name='inflation')
    R_K     ${R^k}$                             (long_name='marginal revenue product of capital')
    q       ${q}$                               (long_name='Tobins q')
    G       ${G}$                               (long_name='Government spending')
    U_C     ${U_C}$                             (long_name='marginal utility wrt consumption')
    U_N     ${U_N}$                             (long_name='marginal disutility of labor')
    mu      ${\mu}$                             (long_name='markup')
    sigma_z     ${\sigma_z}$                    (long_name='TFP shock volatility')
    sigma_G     ${\sigma_G}$                    (long_name='G shock volatility')
    R           ${R}$                           (long_name='policy rate')
    Xi          ${\Xi}$                         (long_name='marginal costs')
    log_Y       ${\hat Y}$      (long_name='log output')
    log_C       ${\hat C}$      (long_name='log consumption')
    log_I       ${\hat I}$      (long_name='log investment')
    log_mu      ${\hat \mu}$    (long_name='log markup')
    log_N       ${\hat N}$      (long_name='log labor')
    log_W       ${\hat W}$      (long_name='log wage')
    log_K       ${\hat K}$      (long_name='log capital stock')
    firm_wedge  ${\tau_f}$      (long_name='firm wedge')
    household_wedge ${\tau_h}$  (long_name='household wedge')
    MPN         ${MP_N}$        (long_name='Marginal product of labor')
    MPK         ${MP_K}$        (long_name='Marginal product of capital')
    MRS         ${MRS}$         (long_name='Marginal rate of substitution')
@#if capital_utilization == 1
    u           ${u}$           (long_name='capital utilization')
@# endif
@#if HP_Taylor == 1
    Y_hp_cyc    ${Y^{hp}_{cyc}}$ (long_name='model consistent HP output gap')
@# endif
    T                   $T$                 (long_name='lump sum taxation')
    pi_annualized       ${\Pi^{ann}}$       (long_name='Annualized Inflation')
    R_annualized        ${R^{ann}}$         (long_name='Annualized gross risk-free interest rate')
    R_K_annualized      ${r^{k,ann}}$       (long_name='Annualized gross return on capital services')
    log_pi_annualized   ${\pi^{ann}}$       (long_name='Log annualized Inflation')
    log_R_annualized    ${R^{ann}}$         (long_name='Log annualized gross risk-free interest rate')
    log_R_K_annualized  ${r^{k,ann}}$       (long_name='Log annualized gross return on capital services')
    Y_net               ${Y^{net}}$         (long_name='Net output')
    M                   ${M}$               (long_name='Stochastic discount factor')
    @#if Separate_Firm_Risk_Aversion
        M_firm          ${M^{firm}}$        (long_name='Stochastic discount factor firm')
    @# endif
;

predetermined_variables K;

varexo eps_z        ${\varepsilon^z}$               (long_name='technology shock')
    eps_G           ${\varepsilon^G}$               (long_name='G shock')
    eps_sigma_z     ${\varepsilon^{\sigma^z}}$      (long_name='technology volatility shock')
    eps_sigma_G     ${\varepsilon^{\sigma^G}}$      (long_name='G volatility shock')
;

parameters siggma   ${\sigma}$      (long_name='risk aversion')
    theta_v         ${\theta_v}$    (long_name='uncertainty resolution preference')
    chi             $\chi$          (long_name='intertemporal elasticity of substitution')
    Frisch_target   $\kappa$        (long_name='Frisch elasticity')
    betta           $\beta$         (long_name='discount factor')
    theta_p         ${\theta_p}$    (long_name='demand elasticity')
    theta_w         ${\theta_w}$    (long_name='wage elasticity')
    labor_share     ${laborshare}$  (long_name='labor share')
    Phi             $\Phi$          (long_name='fixed costs')
    delta           $\delta$        (long_name='depreciation')
    @#if Adjustments_Costs == 1
    phi_k           $\phi_k$        (long_name='capital adjustment costs')
    @# endif
    phi_p           $\phi_P$        (long_name='price adjustment costs')
    phi_w           $\phi_W$        (long_name='wage adjustment costs')
    @#if capital_utilization == 1
        delta_1     ${\delta_1}$    (long_name='linear capital utilization parameter')
        delta_2divdelta_1 ${\frac{\delta_2}{\delta_1}}$ (long_name='capital utilization elasticity')
    @# endif
    Pi_bar          ${\bar \Pi}$    (long_name='steady state inflation')
    rho_r           $\rho_r$        (long_name='interest smoothing')
    phi_Rpi         $\phi_{R\pi}$   (long_name='inflation feedback')
    phi_Ry          $\phi_{Ry}$     (long_name='output growth feedback')
    sigma_z_bar     ${\sigma^z}$        (long_name='Volatility tech shock')
    sigma_G_bar     ${\sigma^G}$        (long_name='Volatility G shock')
    sigma_sigma_z   ${\sigma^{\sigma^z}}$   (long_name='Volatility tech uncertainty shock')
    sigma_sigma_G   ${\sigma^{\sigma^G}}$   (long_name='Volatility G uncertainty shock')
    rho_z           $\rho_z$        (long_name='Autocorrelation tech shock')
    rho_G           $\rho_{G}$      (long_name='Autocorrelation G shock')
    rho_sigma_z     $\rho_{\sigma^z}$       (long_name='Autocorrelation tech uncertainty shock')
    rho_sigma_G     $\rho_{\sigma^{G}}$     (long_name='Autocorrelation G uncertainty shock')
    Z_bar           $Z$             (long_name='SS technology')
    PF_normalization ${Y^{norm}}$   (long_name='PF normalization')
    V_normalization ${V^{norm}}$    (long_name='Value function normalization')
    tau_c           ${\tau^{c}}$    (long_name='Consumption tax')
    tau_l           ${\tau^{l}}$    (long_name='Labor tax')    
    N_overhead      ${N^o}$         (long_name='Overhead labor')
    N_overhead_share $\phi^o$       (long_name='Share of overhead labor')
    gshare          ${\frac{G}{Y}}$ (long_name='Spending to output ratio')
    Gbar            ${\bar G}$      (long_name='Steady state G')
    phi_G_y         ${\phi_{Gy}}$   (long_name='Output feedback on G')
    steady_state_markup_indicator 
    overhead_labor_indicator 
    CES_indicator 
    separate_firm_risk_aversion_indicator
    Iso_elastic_utility_indicator
    @#if Iso_elastic_utility==0
        eta             $\eta$          (long_name='share of consumption in Cobb-Douglas aggregator')
    @#else
        eta             $\eta$          (long_name='labor disutility')
        h               ${h}$           (long_name='habit persistence')
    @# endif
    alppha              ${\alpha}$      (long_name='capital share')
    @#if CES == 1
        psii_CES        ${\psi}$        (long_name='CES substitution elasticity')
    @# endif
    @#if Separate_Firm_Risk_Aversion
        sigma_firm       ${\sigma_{firm}}$      (long_name='Risk aversion firm')
        theta_v_firm     ${\theta_v^{firm}}$    (long_name='uncertainty resolution preference')
        chi_firm         ${\chi_{firm}}$        (long_name='intertemporal elasticity of substitution')
    @# endif
;

tau_c=0.094;
tau_l=0.22;

betta =0.995;
delta =0.025;
gshare=0.205;

@#if Adjustments_Costs 
    @#if CEE == 0
        phi_k = 2.09;
    @# else
        %CEE
        @#if FV_et_al_calibration==0
            phi_k = 2.5;     
        @# else
            phi_k = 0.75;
        @# endif        
    @# endif
@# endif

@#if capital_utilization == 1
    delta_1=0; %set in steady state
    delta_2divdelta_1=0.02/(1/betta-(1-delta));
@# endif

@#if overhead_labor == 1
    N_overhead_share=0.11;
    N_overhead=0; %set in steady state
    overhead_labor_indicator=1;
@# else
    N_overhead_share=0;
    N_overhead=0;
    overhead_labor_indicator=0;
@# endif

@#if CES == 1
    psii_CES=0.5;
    labor_share=2/3;
    alppha =1/3; %set in steady state file based on N_overhead and labor_share
    CES_indicator=1;
@# else
    labor_share=2/3;
    alppha =1/3; %set in steady state file based on N_overhead and labor_share
    CES_indicator=0;
@# endif

@#if FV_et_al_calibration==0
    theta_p = 11; %demand elasticity
    Calvo_duration=2/3;
@# else
    theta_p = 21; %demand elasticity
    Calvo_duration=0.75;
@# endif        
@#if sticky_prices == 1
    phi_p =(theta_p-1)*Calvo_duration/((1-Calvo_duration)*(1-betta*Calvo_duration)); %0.75 corresponds to 4 quarter duration
@# else
    phi_p =0;
@# endif

@#if FV_et_al_calibration==0
    theta_w = 11; %demand elasticity
    Calvo_wage_target=2/3;
@# else
    theta_w = 21; %demand elasticity
    Calvo_wage_target=0.75;
@# endif        
Frisch_target=1;
@#if sticky_wages == 1
    phi_w = (theta_w-1)*(1-tau_l)*labor_share*Calvo_wage_target*(1+theta_w*Frisch_target)/((1-Calvo_wage_target)*(1-betta*Calvo_wage_target));  %0.75 corresponds to 4 quarter duration
    steady_state_markup_indicator=1; %indicator used in steady state file
@# else
    phi_w =0;
    steady_state_markup_indicator=1;
@# endif


Pi_bar =1;
rho_r =0.75;
phi_Rpi =1.35;
phi_Ry =0.125;


chi = 0.5;
separate_firm_risk_aversion_indicator=0;
@#if CRRA == 1
    siggma = 1/chi;
    @#if Separate_Firm_Risk_Aversion
        sigma_firm=20;
        chi_firm=1/sigma_firm;
        separate_firm_risk_aversion_indicator=1;
    @# endif
@# else
    %EZ
    siggma = 66;
    @#if Separate_Firm_Risk_Aversion
        sigma_firm=66;
        chi_firm=1/sigma_firm;
        separate_firm_risk_aversion_indicator=1;
    @# endif
@# endif

@#if Iso_elastic_utility==0
    Iso_elastic_utility_indicator=0;
@# else
    Iso_elastic_utility_indicator=1;
    eta=0; %set in SS file
    h=0.7;
@# endif

Z_bar =1;

%% Mean parameters from simulation_2016_12_10_10_59
%% No exp()-formulation with HP-trend and contemporaneous effect of volatility

@#if Leduc_Liu_calibration==0
    rho_z=0.773;
    rho_sigma_z=0.517;
    sigma_z_bar =0.007;
    sigma_sigma_z= 0.0023;
@#else
    rho_z=0.95;
    rho_sigma_z=0.76;
    sigma_z_bar =0.01;
    sigma_sigma_z= 0.005;
@#endif


%% Mean parameters from simulation_2020_2_16_15_27
%% No exp()-formulation with output feedback
phi_G_y=0.022;
rho_G =0.938;
rho_sigma_G =0.504 ;
sigma_G_bar =0.008;
sigma_sigma_G =0.003;

eta=0; %set in SS file 

Phi = 0; %fixed costs, set in SS file

model;
#log_R_bar=log(Pi_bar/betta);
@#if capital_utilization == 1
    #delta_2=delta_2divdelta_1*delta_1;
@# endif

@#if Iso_elastic_utility==0
    [name='Definition value function']
    V=(V_normalization*(C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v) +betta*(E_t_V_tp1_1_minus_sigma)^(1/theta_v))^(theta_v/(1-siggma));
    [name='Definition of marginal utility w.r.t. to consumption']
    U_C=V^(1-(1-siggma)/theta_v)*eta*V_normalization*(C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v)/C;

    [name='Definition of marginal utility w.r.t. to hours worked']
    U_N=-(V_normalization*(C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v)+betta*(E_t_V_tp1_1_minus_sigma)^(1/theta_v))^(theta_v/(1-siggma)-1)*V_normalization*(C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v)*(1-eta)/(1-N);
    
    @#if Separate_Firm_Risk_Aversion
        [name='Definition SDF firms']
        M_firm=betta*((C^eta*(1-N)^(1-eta))/(C(-1)^eta*(1-N(-1))^(1-eta)))^((1-sigma_firm)/theta_v_firm)*(C(-1)/C)*(V^(1-sigma_firm)/E_t_V_tp1_1_minus_sigma(-1))^(1-1/theta_v_firm);
    @# endif

    [name='Definition SDF']
    M=betta*((C^eta*(1-N)^(1-eta))/(C(-1)^eta*(1-N(-1))^(1-eta)))^((1-siggma)/theta_v)*(C(-1)/C)*(V^(1-siggma)/E_t_V_tp1_1_minus_sigma(-1))^(1-1/theta_v);
@#else
    [name='Definition value function']
    V=V_normalization*((C-h*C(-1))^(1-siggma)/(1-siggma)-eta*N^(1+Frisch_target)/(1+Frisch_target)) +betta*V(+1);

    [name='Definition of marginal utility w.r.t. to consumption']
    U_C=V_normalization*(C-h*C(-1))^(-siggma);

    [name='Definition of marginal utility w.r.t. to hours worked']
    U_N=-V_normalization*eta*N^(Frisch_target);
    
    @#if Separate_Firm_Risk_Aversion
        [name='Definition SDF firms']
        M_firm=betta*((C-h*C(-1))^(-sigma_firm))/((C(-1)^-h*C(-2))^(-sigma_firm));
    @# endif

    [name='Definition SDF']
    M=betta*((C-h*C(-1))^(-siggma))/((C(-1)-h*C(-2))^(-siggma));
@#endif

[name='Auxiliary variable needed for Epstein-Zin preferences']
E_t_V_tp1_1_minus_sigma=V(+1)^(1-siggma);

[name='Budget constraint']
(1+tau_c)*C=(1-tau_l)*W*N+D-phi_w/2*(W/W(-1)*Pi/Pi_bar-1)^2*Y+T;

[name='labor FOC']
U_N*(-theta_w)*N + 
U_C/(1+tau_c)*((1-theta_w)*(1-tau_l)*W*N-phi_w*(W/W(-1)*Pi/Pi_bar-1)*W/W(-1)*Pi/Pi_bar*Y)
+betta*U_C(+1)/(1+tau_c(+1))*(phi_w*(W(+1)/W*Pi(+1)/Pi_bar-1)*W(+1)/W*Pi(+1)/Pi_bar*Y(+1))=0;        


[name='Production function']

@#if CES == 1
    @#if capital_utilization == 1
        Y=PF_normalization*(alppha*(K*u)^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES)-Phi;
    @#else
        Y=PF_normalization*(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES)-Phi;
    @# endif
@# else
    @#if capital_utilization == 1
        Y=PF_normalization*(K*u)^alppha*(Z*(N-N_overhead))^(1-alppha)-Phi;
    @#else
        Y=PF_normalization*K^alppha*(Z*(N-N_overhead))^(1-alppha)-Phi;
    @# endif
@# endif

@#if Adjustments_Costs
    @#if CEE == 0
    %Hayashi
        @#if capital_utilization == 1
            [name='LOM Capital']
            K(+1)=(1-(delta+delta_1*(u-1)+delta_2/2*(u-1)^2)-phi_k/2*(I/K-delta)^2)*K+I;
            [name='FOC capital']
            q=M(+1)*(R_K(+1)*u(+1)+q(+1)*(1-(delta+delta_1*(u(+1)-1)+delta_2/2*(u(+1)-1)^2)-phi_k/2*(I(+1)/K(+1)-delta)^2 +phi_k*(I(+1)/K(+1)-delta)*(I(+1)/K(+1))));
            [name='FOC investment']
            1/q=1-phi_k*(I/K-delta);
            [name='FOC utilization']
            R_K=q*(delta_1+delta_2*(u-1));
        @#else
            [name='LOM Capital']
            K(+1)=(1-delta-phi_k/2*(I/K-delta)^2)*K+I;
            [name='FOC capital']
            q=M(+1)*(R_K(+1)+q(+1)*(1-delta-phi_k/2*(I(+1)/K(+1)-delta)^2 +phi_k*(I(+1)/K(+1)-delta)*(I(+1)/K(+1))));
            [name='FOC investment']
            1/q=1-phi_k*(I/K-delta);
        @# endif
    @# else
    %CEE
        @#if capital_utilization == 1
            [name='LOM Capital']
            K(+1)=(1-(delta+delta_1*(u-1)+delta_2/2*(u-1)^2))*K+I*(1-phi_k/2*(I/I(-1)-1)^2);
            [name='FOC capital']
            q = M(+1)*(R_K(+1)*u(+1)+q(+1)*(1-(delta+delta_1*(u(+1)-1)+delta_2/2*(u(+1)-1)^2)));
            [name='FOC investment']
            1 = q*(1-phi_k/2*(I/I(-1)-1)^2-phi_k*(I/I(-1)-1)*I/I(-1))
                +M(+1)*q(+1)*phi_k*(I(+1)/I-1)*(I(+1)/I)^2;
            [name='FOC utilization']            
            R_K= q*(delta_1+delta_2*(u-1));
        @#else
            [name='LOM Capital']
            K(+1)=(1-delta)*K+I*(1-phi_k/2*(I/I(-1)-1)^2);
            [name='FOC capital']
            q = M(+1)*(R_K(+1)+q(+1)*(1-delta));
            [name='FOC investment']
            1 = q*(1-phi_k/2*(I/I(-1)-1)^2-phi_k*(I/I(-1)-1)*I/I(-1))
                +M(+1)*q(+1)*phi_k*(I(+1)/I-1)*(I(+1)/I)^2;
        @# endif
    @# endif
@# else
    [name='LOM Capital']
    K(+1)=(1-delta)*K+I;
    [name='FOC capital']
    q = M(+1)*(R_K(+1)+q(+1)*(1-delta));
    [name='FOC investment']
    1 = q;
@# endif


[name='Cash flows']
D=Y-W*N-I-phi_p/2*(Pi/Pi_bar-1)^2*Y;

[name='Labor FOC']
W=Xi*MPN;

[name='Capital FOC']
@#if capital_utilization == 1
    R_K=Xi*MPK/u;
@#else
    R_K=Xi*MPK;    
@# endif
   
@#if CES == 1
    @#if capital_utilization == 1
        [name='Marginal product of labor']
        MPN=PF_normalization*(alppha*(K*u)^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES-1)*(1-alppha)*(Z*(N-N_overhead))^psii_CES/(N-N_overhead);
        [name='Marginal product of caital']
        MPK=PF_normalization*(alppha*(K*u)^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES-1)*alppha*(u)^psii_CES*K^(psii_CES-1);
    @#else
        [name='Marginal product of labor']
        MPN=PF_normalization*(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES-1)*(1-alppha)*(Z*(N-N_overhead))^psii_CES/(N-N_overhead);
        [name='Marginal product of caital']
        MPK=PF_normalization*(alppha*K^psii_CES+(1-alppha)*(Z*(N-N_overhead))^psii_CES)^(1/psii_CES-1)*alppha*K^(psii_CES-1);
    @# endif
@# else
    @#if capital_utilization == 1
        [name='Marginal product of labor']
        MPN=(1-alppha)*PF_normalization*(K*u)^alppha*(Z*(N-N_overhead))^(-alppha)*Z;
        [name='Marginal product of caital']
        MPK=alppha*PF_normalization*(K*u)^(alppha)/K*(Z*(N-N_overhead))^(1-alppha);
    @#else
        [name='Marginal product of labor']
        MPN=(1-alppha)*PF_normalization*K^alppha*(Z*(N-N_overhead))^(-alppha)*Z;
        [name='Marginal product of caital']
        MPK=alppha*PF_normalization*K^(alppha-1)*(Z*(N-N_overhead))^(1-alppha);
    @# endif
@# endif

@#if Separate_Firm_Risk_Aversion
    @#if FVEtal_Rotemberg
        [name='Pricing FOC']
        phi_p*(Pi/Pi_bar-1)*(Pi/Pi_bar)=theta_p*phi_p/2*(Pi/Pi_bar-1)^2+(1-theta_p)+theta_p*Xi+
            phi_p*M_firm(+1)*(Y(+1)/Y)*(Pi(+1)/Pi_bar-1)*(Pi(+1)/Pi_bar);           
    @# else
        [name='Pricing FOC']
        phi_p*(Pi/Pi_bar-1)*(Pi/Pi_bar)=(1-theta_p)+theta_p*Xi+
            phi_p*M_firm(+1)*(Y(+1)/Y)*(Pi(+1)/Pi_bar-1)*(Pi(+1)/Pi_bar);           
    @# endif        
@# else
    @#if FVEtal_Rotemberg
        [name='Pricing FOC']
        phi_p*(Pi/Pi_bar-1)*(Pi/Pi_bar)=theta_p*phi_p/2*(Pi/Pi_bar-1)^2+(1-theta_p)+theta_p*Xi+
            phi_p*M(+1)*(Y(+1)/Y)*(Pi(+1)/Pi_bar-1)*(Pi(+1)/Pi_bar);           
    @# else
        [name='Pricing FOC']
        phi_p*(Pi/Pi_bar-1)*(Pi/Pi_bar)=(1-theta_p)+theta_p*Xi+
            phi_p*M(+1)*(Y(+1)/Y)*(Pi(+1)/Pi_bar-1)*(Pi(+1)/Pi_bar);           
    @# endif        
@# endif

[name='Bond FOC']
1=M(+1)*(R/Pi(+1));

[name='Markup definition']
mu=1/Xi;

[name='Gross output']
Y_net=(1-phi_p/2*(Pi/Pi_bar-1)^2-phi_w/2*(W/W(-1)*Pi/Pi_bar-1)^2)*Y;

[name='Definition marginal rate of substitution']
MRS=-U_N/U_C; %equal to (1-eta)/eta*C/(1-N)

[name='Definition firm wedge']
firm_wedge=log(MPN/W);

[name='Definition household wedge']
household_wedge=-log(MRS*((1+tau_c)/(1-tau_l))/W); 
 
[name='Taylor Rule']
@#if HP_Taylor == 0
    log(R)=rho_r*log(R(-1))+(1-rho_r)*(log_R_bar+phi_Rpi*(log(Pi)-log(Pi_bar))+phi_Ry*(log(Y)-log(Y(-1))))
;
@# else
    log(R)=rho_r*log(R(-1))+(1-rho_r)*(log_R_bar+phi_Rpi*(log(Pi)-log(Pi_bar))+phi_Ry*Y_hp_cyc);
    [name='HP-filtered output gap']
    Y_hp_cyc*(1+6*1600)+Y_hp_cyc(-1)*(-4*1600)+Y_hp_cyc(+1)*(-4*1600)+Y_hp_cyc(-2)*1600+Y_hp_cyc(+2)*1600=Y*(6*1600)+Y(-1)*(-4*1600)+Y(+1)*(-4*1600)+Y(-2)*1600+Y(+2)*1600;
@# endif

[name='Government budget constraint']
tau_c*C+tau_l*W*N=T+Gbar*G;

%Shock processes
Z=(1-rho_z)*Z_bar+rho_z*Z(-1)+sigma_z*eps_z;
sigma_z=(1-rho_sigma_z)*sigma_z_bar+rho_sigma_z*sigma_z(-1)+sigma_sigma_z*eps_sigma_z;
G=(1-rho_G)*1+rho_G*G(-1)+phi_G_y*(log(Y(-1)/steady_state(Y)))+sigma_G*eps_G;
sigma_G=(1-rho_sigma_G)*sigma_G_bar+rho_sigma_G*sigma_G(-1)+sigma_sigma_G*eps_sigma_G;

%Observed variables
log_Y=log(Y_net);
log_C=log(C);
log_I=log(I);
log_mu=log(mu);
log_N=log(N);
log_K=log(K);
log_W =log(W);
pi_annualized=1+(4*(Pi-1));
R_annualized=1+(4*(R-1));
R_K_annualized=4*R_K;
log_pi_annualized=log(pi_annualized);
log_R_annualized=log(R_annualized);
log_R_K_annualized=log(R_K_annualized);

end;


shocks;
var eps_z; stderr 1;
var eps_G; stderr 1;
var eps_sigma_z ; stderr 1;
var eps_sigma_G; stderr 1;
end;

% steady;

% write_latex_original_model;
% write_latex_static_model;
% write_latex_definitions;
% write_latex_parameter_table;
% collect_latex_files;
%       if system(['pdflatex -halt-on-error -interaction=batchmode ' M_.fname '_TeX_binder.tex'])
%           error('TeX-File did not compile.')
%       end
% check;

verbatim;
shocks={'eps_sigma_z',
        'eps_sigma_G'
        'eps_z'
        'eps_G'};
@#if Large_shocks == 0
shock_size=[2 0 0 0;
            0 2 0 0 
            0 0 1 0
            0 0 0 1];
@#else
shock_size=[4 0 0 0;
            0 4 0 0 
            0 0 1 0
            0 0 0 1];
@#endif
end;

stoch_simul(order=4,k_order_solver,nograph,noprint,nomoments,nofunctions,nocorr,irf=0) log_Y log_K log_C log_I log_mu log_N sigma_z sigma_G log_pi_annualized log_R_annualized log_W firm_wedge household_wedge G Z;

burnin=500;
options_.irf=20;

for shock_iter=1:size(shocks,1);
    ex_=zeros(burnin+options_.irf,M_.exo_nbr);
    @#if Sequential_shocks == 1
        ex_(burnin+1,strmatch(shocks{shock_iter},M_.exo_names,'exact'))=shock_size(shock_iter,shock_iter)/2;
        ex_(burnin+2,strmatch(shocks{shock_iter},M_.exo_names,'exact'))=shock_size(shock_iter,shock_iter)/2;
    @#else
        ex_(burnin+1,strmatch(shocks{shock_iter},M_.exo_names,'exact'))=shock_size(shock_iter,shock_iter);
    @#endif
    IRF_mat=simult_(M_,options_,oo_.dr.ys,oo_.dr,ex_,options_.order)';

    stochastic_steady_state=IRF_mat(burnin,:); % stochastic_steady_state/EMAS is any of the final points after burnin


    IRF_mat_percent_from_SSS_non_logged = (IRF_mat(1+burnin+1:1+burnin+options_.irf,:)-repmat(stochastic_steady_state,options_.irf,1))./repmat(stochastic_steady_state,options_.irf,1); %only valid for variables not yet logged
    IRF_mat_percent_from_SSS_logged = (IRF_mat(1+burnin+1:1+burnin+options_.irf,:)-stochastic_steady_state);


    var_names={'sigma_z','sigma_G'};
    for var_iter=1:length(var_names)
        var_pos 	= strmatch(var_names{var_iter},M_.endo_names,'exact');
        oo_.irfs.([var_names{var_iter},'_',shocks{shock_iter}])=IRF_mat_percent_from_SSS_non_logged(:,var_pos);
    end
    var_names={'firm_wedge','household_wedge','Z','G','log_Y','log_C','log_I','log_N','log_W','log_K',... 
    'log_pi_annualized','log_R_annualized'};

    for var_iter=1:length(var_names)
        var_pos 	= strmatch(var_names{var_iter},M_.endo_names,'exact');
        oo_.irfs.([var_names{var_iter},'_',shocks{shock_iter}])=IRF_mat_percent_from_SSS_logged(:,var_pos);
    end
end
                
Create_1_by_4_vola_figures

var_names={'sigma_G','log_Y','log_C','log_I','log_N','log_W','log_K',...'firm_wedge','household_wedge',...
    'log_pi_annualized','log_R_annualized'};

h_sigma_G=figure('Name','IRFs Government Spending Volatility Shock');

set(h_sigma_G,'PaperType','a4','PaperPositionMode','manual','PaperUnits','centimeters','PaperPosition',[1,1,28,21],'renderer', 'painters');
for var_iter=1:length(var_names)
    subplot(3,3,var_iter)
    if 100*max(abs(oo_.irfs.([var_names{var_iter},'_eps_sigma_G']))) >1e-5
        plot(1:options_.irf,100*oo_.irfs.([var_names{var_iter},'_eps_sigma_G']),'LineWidth',1.5)
    else
        plot(1:options_.irf,zeros(size(oo_.irfs.([var_names{var_iter},'_eps_sigma_G']))),'LineWidth',1.5)
    end
    temp=deblank(M_.endo_names_tex{strmatch(var_names{var_iter},M_.endo_names,'exact')});
    title(['$' temp(2:end-1) '$'] ,'FontSize',10,'Interpreter','Latex')
    axis tight
    if var_iter>7
        ylabel('p.p.','FontSize',10)
    else
        ylabel('percent','FontSize',10)
    end
    if var_iter>6
        xlabel('quarters','FontSize',10)        
    end
    set(gca,'FontSize',10)
    hline(0);
end

var_names={'G','log_Y','log_C','log_I','log_N','log_W','log_K',...'firm_wedge','household_wedge',...
    'log_pi_annualized','log_R_annualized'};

h_G=figure('Name','IRFs Government Spending Shock');

set(h_G,'PaperType','a4','PaperPositionMode','manual','PaperUnits','centimeters','PaperPosition',[1,1,28,21],'renderer', 'painters');
for var_iter=1:length(var_names)
    subplot(3,3,var_iter)
    if 100*max(abs(oo_.irfs.([var_names{var_iter},'_eps_G']))) >1e-5
        plot(1:options_.irf,100*oo_.irfs.([var_names{var_iter},'_eps_G']),'LineWidth',1.5)
    else
        plot(1:options_.irf,zeros(size(oo_.irfs.([var_names{var_iter},'_eps_G']))),'LineWidth',1.5)
    end
    temp=deblank(M_.endo_names_tex{strmatch(var_names{var_iter},M_.endo_names,'exact')});
    title(['$' temp(2:end-1) '$'] ,'FontSize',10,'Interpreter','Latex')
    axis tight
    if var_iter>7
        ylabel('p.p.','FontSize',10)
    else
        ylabel('percent','FontSize',10)
    end
    if var_iter>6
        xlabel('quarters','FontSize',10)        
    end
    set(gca,'FontSize',10)
    hline(0);
end

for var_iter=2:length(var_names)
    var_pos 	= strmatch(var_names{var_iter},M_.endo_names,'exact');
    oo_.irfs.([var_names{var_iter},'_',shocks{shock_iter}])=IRF_mat_percent_from_SSS_logged(:,var_pos);
end

h_sigma_z=figure('Name','IRFs Technology Volatility Shock');

set(h_sigma_z,'PaperType','a4','PaperPositionMode','manual','PaperUnits','centimeters','PaperPosition',[1,1,28,21],'renderer', 'painters');
for var_iter=1:length(var_names)
    subplot(3,3,var_iter)
    if 100*max(abs(oo_.irfs.([var_names{var_iter},'_eps_sigma_z']))) >1e-5
        plot(1:options_.irf,100*oo_.irfs.([var_names{var_iter},'_eps_sigma_z']),'LineWidth',1.5)
    else
        plot(1:options_.irf,zeros(size(oo_.irfs.([var_names{var_iter},'_eps_sigma_z']))),'LineWidth',1.5)
    end
    temp=deblank(M_.endo_names_tex{strmatch(var_names{var_iter},M_.endo_names,'exact')});
    title(['$' temp(2:end-1) '$'] ,'FontSize',10,'Interpreter','Latex')
    axis tight
    if var_iter>7
        ylabel('p.p.','FontSize',10)
    else
        ylabel('percent','FontSize',10)
    end
    if var_iter>6
        xlabel('quarters','FontSize',10)        
    end
    set(gca,'FontSize',10)
    hline(0);
end


var_names={'Z','log_Y','log_C','log_I','log_N','log_W','log_K',... 'firm_wedge','household_wedge',...
    'log_pi_annualized','log_R_annualized'};

h_Z=figure('Name','IRFs Technology Shock');

set(h_Z,'PaperType','a4','PaperPositionMode','manual','PaperUnits','centimeters','PaperPosition',[1,1,28,21],'renderer', 'painters');
for var_iter=1:length(var_names)
    subplot(3,3,var_iter)
    if 100*max(abs(oo_.irfs.([var_names{var_iter},'_eps_z']))) >1e-5
        plot(1:options_.irf,100*oo_.irfs.([var_names{var_iter},'_eps_z']),'LineWidth',1.5)
    else
        plot(1:options_.irf,zeros(size(oo_.irfs.([var_names{var_iter},'_eps_z']))),'LineWidth',1.5)
    end
    temp=deblank(M_.endo_names_tex{strmatch(var_names{var_iter},M_.endo_names,'exact')});
    title(['$' temp(2:end-1) '$'] ,'FontSize',10,'Interpreter','Latex')
    axis tight
    if var_iter>7
        ylabel('p.p.','FontSize',10)
    else
        ylabel('percent','FontSize',10)
    end
    if var_iter>6
        xlabel('quarters','FontSize',10)        
    end
    set(gca,'FontSize',10)
    hline(0);
end
