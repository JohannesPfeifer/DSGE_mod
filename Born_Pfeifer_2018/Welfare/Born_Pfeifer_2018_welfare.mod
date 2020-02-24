/*
 * This file implements the mod-file for the Welfare analysis in Born/Pfeifer (2018):
 * "The New Keynesian Wage Phillips Curve: Calvo vs. Rotemberg". 
 * It contains the baseline New Keynesian model of Galí (2015), Chapter 6 in nonlinear form
 * with the EHL and SGU insurance schemes and Calvo and Rotemberg wage setting, while price
 * setting follows Calvo.
 *
 * It computes consumption equivalents relative to the flex-price natural allocation as described
 * in the paper. For doing so, it uses a solver to find the correct consumption equivalent lambda
 * based on the recursive lifetime welfare measures defined within the model block.
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model. 
 */

/*
 * Copyright (C) 2018 Johannes Pfeifer and Benjamin Born
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
@#ifndef strict_targeting
    @#define strict_targeting=1
@#endif        

@#ifndef price_targeting
    @#define price_targeting=1
@#endif        
@#ifndef wage_targeting
    @#define wage_targeting=0
@#endif        
@#ifndef composite_targeting
    @#define composite_targeting=0
@#endif           
        
@#ifndef Ramsey_policy
    @#define Ramsey_policy=0
@#endif        


@#ifndef sticky_prices
    @#define sticky_prices=1
@#endif        

@#ifndef sticky_wages
    @#define sticky_wages=1
@#endif        

@#if Ramsey_policy==0
    @#if strict_targeting == 1
        table_title='Strict targeting';
    @#else
        table_title='Flexible targeting';
    @#endif

    @#if price_targeting ==1
         case_title='Price';
    @#else
        @#if wage_targeting ==1
            case_title='Wage';
        @#else
            @#if composite_targeting ==1
                case_title='Composite';
            @#else
                error('One case must be set to 1')                
            @#endif
        @#endif       
    @#endif
@#else
    table_title='Ramsey Policy';
    case_title='Optimal';
@#endif

@#ifndef efficient_steady_state
    @#define efficient_steady_state=1
@#endif        

@#ifndef log_utility
    @#define log_utility=1
@#endif        
    
@#ifndef fixed_WPC_slope
    %compute parameters implied by WPC slope instead of taking Calvo and
    %computing implied Rotemberg
    @#define fixed_WPC_slope=0
@#endif        

@#ifndef SGU_framework
    @#define SGU_framework=0
@#endif        

@#ifndef Calvo
    @#define Calvo=1
@#endif        
@#if Calvo == 1
    table_title=[table_title,' Calvo'];
@#else
    table_title=[table_title,' Rotemberg'];
@#endif


@#ifndef taxes
    @#define taxes=0
@#endif        

@#ifndef compare_linear_model
    @#define compare_linear_model=0
@#endif        

var C               ${C}$           (long_name='Consumption')
    W_real          ${\frac{W}{P}}$ (long_name='Real Wage')
    Pi              ${\Pi}$         (long_name='inflation')
    A               ${A}$           (long_name='AR(1) technology process')
    N               ${N}$           (long_name='Hours worked')
    R               ${R^n}$         (long_name='Nominal Interest Rate') 
    realinterest    ${R^{r}}$       (long_name='Real Interest Rate')
    Y               ${Y}$           (long_name='Output') 
    Q               ${Q}$           (long_name='Bond price')
    Z               ${Z}$           (long_name='AR(1) preference shock process')
    S               ${S}$           (long_name='Price dispersion')
    @#if sticky_prices == 1
        Pi_star         ${\Pi^*}$       (long_name='Optimal reset price')
        x_aux_1         ${x_1}$         (long_name='aux. var. 1 recursive price setting')
        x_aux_2         ${x_2}$         (long_name='aux. var. 2 recursive price setting')
    @#endif   
    MC              ${mc}$          (long_name='real marginal costs')
    M_real          ${M/P}$         (long_name='real money stock')
    i_ann           ${i^{ann}}$     (long_name='annualized nominal interest rate')
    pi_p_ann        ${\pi^{ann}}$   (long_name='annualized inflation rate')
    r_real_ann      ${r^{r,ann}}$   (long_name='annualized real interest rate')
    log_y           ${\log(Y)}$      (long_name='log output')
    log_W_real      ${\log(W/P)}$    (long_name='log real wage')
    log_N           ${\log(N)}$      (long_name='log hours')
    log_A           ${\log(A)}$      (long_name='log technology level')
    log_Z           ${\log(Z)}$      (long_name='log preference shock')
    log_Pi          ${\log(\Pi)}$    (long_name='log inflation')
    log_Pi_w        ${\log(\Pi_w)}$  (long_name='log wage inflation')
    nu              ${\nu}$         (long_name='AR(1) monetary policy shock process')    
    Pi_w            ${\Pi^w}$       (long_name='wage inflation')
    @#if Calvo
        Pi_w_star       ${\Pi^*_w}$ (long_name='Optimal reset wage')
        f               ${f}$       (long_name='aux. var. recursive wage setting')
        S_w             ${S^w}$     (long_name='Wage dispersion')
    @#endif
    pi_w_ann        ${\pi^{w,ann}}$ (long_name='annualized wage inflation rate')
    N_d             ${N^d}$         (long_name='Effective Hours worked')
    tau_n           ${\tau_n}$      (long_name='labor tax')
    MPN             ${MPN}$         (long_name='Marginal product of labor')
    MRS             ${MRS}$         (long_name='After tax marginal rate of substitution')
    Utility         ${\mathbb{W}}$  (long_name='Period Utility')
    X_w             ${X^{W}}$       (long_name='Auxiliary variable from aggregation in utility')
    N_nat           ${N^{nat}}$     (long_name='Natural hours worked')
    Y_nat           ${Y^{nat}}$     (long_name='Natural output')
    log_Y_gap       ${\log{Y^{gap}}}$  (long_name='Log output gap')
    @#if Ramsey_policy==0    
        Recursive_Welfare  ${\mathbb{V}}$  (long_name='Recursive welfare')
        Welfare_gap        ${\mathbb{V^{gap}}}$  (long_name='Welfare gap given lambda')
    @#endif    
    Recursive_natural_welfare_equivalent ${\mathbb{V^{nat}_{gap}}}$  (long_name='Recursive welfare natural EQ, given lambda')
    ;

varexo eps_a        ${\varepsilon_a}$   (long_name='technology shock')
       eps_z        ${\varepsilon_z}$   (long_name='preference shock')
       eps_nu       ${\varepsilon_\nu}$ (long_name='monetary policy shock')
       eps_tau_n    ${\varepsilon_{\tau^n}}$ (long_name='labor tax shock')
       ;   

parameters alppha       ${\alpha}$      (long_name='capital share')
    betta               ${\beta}$       (long_name='discount factor')
    rho_a               ${\rho_a}$      (long_name='autocorrelation technology shock')
    rho_nu              ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')
    rho_z               ${\rho_{z}}$    (long_name='autocorrelation demand shock')
    siggma              ${\sigma}$      (long_name='inverse EIS')
    varphi              ${\varphi}$     (long_name='inverse Frisch elasticity')
    phi_pi              ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    phi_y               ${\phi_{y}}$    (long_name='output feedback Taylor Rule')
    eta                 ${\eta}$        (long_name='semi-elasticity of money demand')
    epsilon_p           ${\epsilon_p}$  (long_name='demand elasticity')
    theta_p             ${\theta_p}$    (long_name='Calvo parameter')
    tau_s_p             ${\tau^{s}_p}$  (long_name='labor subsidy firms')
    tau_s_w             ${\tau^{s}_w}$  (long_name='labor subsidy workers')
    epsilon_w           ${\epsilon_w}$  (long_name='demand elasticity labor services')
    theta_w             ${\theta_w}$    (long_name='Calvo parameter wages')
    phi_w               ${\phi_w}$      (long_name='Rotemberg parameter wages')
    tau_n_SS            ${\bar {\tau^n}}$   (long_name='Steady state labor tax')
    tau_c               ${\tau_c}$          (long_name='Consumption taxes')
    rho_tau_n           ${\rho_{\tau^n}}$   (long_name='labor tax persistence')
    lambda_w            ${\lambda_{w}}$     (long_name='slope of the Wage PC')
    lambda_p            ${\lambda_{p}}$     (long_name='slope of the PC')
    lambda_utility      ${\lambda}$         (long_name='Welfare_relevant consumption fraction')
    ;
    
%----------------------------------------------------------------
% Parametrization, p. 67  and p. 113-115
%----------------------------------------------------------------
@#if log_utility
    siggma = 1;
@#else
    siggma = 2;
@#endif

varphi=5;
phi_pi = 1.5;
phi_y  = 0.125;
theta_p=0.75;
rho_nu =0.5;
rho_z  = 0.5;
rho_a  = 0.9;
betta  = 0.99;
eta  =3.77; %footnote 11, p. 115
alppha=1/4;
epsilon_p=9;
epsilon_w=4.5;
@#if efficient_steady_state
    tau_s_p=(epsilon_p-1)/epsilon_p;
    tau_s_w=epsilon_w/(epsilon_w-1);
@#else
    tau_s_p=1;
    tau_s_w=1;
@#endif


@#if taxes
    tau_c=0.07;
    tau_n_SS=0.2;
    rho_tau_n=0.9;
@#else
    tau_c=0;
    tau_n_SS=0;
    rho_tau_n=0;
@#endif

@#if fixed_WPC_slope
    %compute parameters implied by WPC slope instead of taking Calvo and
    %computing implied Rotemberg
    lambda_w=0.003652482269504; %0.75 Calvo in EHL
    theta_w=0; %not used with Rotemberg, but should be initialized in any case 
    phi_w=0; %not used with Calvo, but should be initialized in any case 
@#else
    theta_w=3/4;
@#endif     

lambda_utility=0; % fraction of natural consumption required to make agent indifferent (used in optimizer below)
%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model;
@#if Calvo
    [name='LOM wages']
    1=theta_w*Pi_w^(epsilon_w-1)+(1-theta_w)*(Pi_w_star)^(1-epsilon_w);
    [name='LOM wage dispersion']
    S_w=(1-theta_w)*Pi_w_star^(-epsilon_w)+theta_w*Pi_w^(epsilon_w)*S_w(-1);
    @#if SGU_framework
        [name='Auxiliary wage setting recursion 1']
        f=(epsilon_w-1)/epsilon_w*tau_s_w*(1-tau_n)/(1+tau_c)*Z*C^(-siggma)*W_real*Pi_w_star^(1-epsilon_w)*N_d
                +betta*theta_w*(Pi_w_star(+1)/Pi_w_star*Pi_w(+1))^(epsilon_w-1)*f(+1);
        [name='Auxiliary wage setting recursion 2']
        f=Z*Pi_w_star^(-epsilon_w)*N_d^(1+varphi)
                +betta*theta_w*(Pi_w_star(+1)/Pi_w_star*Pi_w(+1))^(epsilon_w)*f(+1);
    @#else
        [name='Auxiliary wage setting recursion 1']
        f=(epsilon_w-1)/epsilon_w*tau_s_w*(1-tau_n)/(1+tau_c)*Z*C^(-siggma)*W_real*Pi_w_star^(1-epsilon_w)*N_d
                +betta*theta_w*(Pi_w_star(+1)/Pi_w_star*Pi_w(+1))^(epsilon_w-1)*f(+1);
        [name='Auxiliary wage setting recursion 2']
        f=Z*Pi_w_star^(-epsilon_w*(1+varphi))*N_d^(1+varphi)
                +betta*theta_w*(Pi_w_star(+1)/Pi_w_star*Pi_w(+1))^(epsilon_w*(1+varphi))*f(+1);
    @#endif
    [name='Definition effective labor']
    N_d=N/S_w;    
@#else
    [name='FOC wage setting']
    0=epsilon_w*(1+tau_c)*N^varphi/C^(-siggma)/W_real+(1-epsilon_w)*tau_s_w*(1-tau_n)-phi_w*(Pi_w/steady_state(Pi)-1)*Pi/steady_state(Pi)*Y/(W_real*N)
        +betta*Z(+1)*C(+1)^(-siggma)/(Z*C^(-siggma))/(W_real*N)*(phi_w*(Pi_w(+1)/steady_state(Pi)-1)*Pi_w(+1)/steady_state(Pi)*Y(+1));
    [name='Definition effective labor']
    N_d=N;    
@#endif
    [name='Definition wage inflation']    
    Pi_w=W_real/W_real(-1)*Pi;            
    [name='Euler equation eq. (3)']
    Q=betta*(C(+1)/C)^(-siggma)*(Z(+1)/Z)/Pi(+1);
    [name='Definition nominal interest rate), p. 22 top']
    R=1/Q;
    [name='Aggregate output, above eq. (14)']
    Y=A*(N_d/S)^(1-alppha);
    [name='Definition Real interest rate']
    R=realinterest*Pi(+1);
    @#if Ramsey_policy==0
        [name='Simple rule']
        @#if strict_targeting == 1
            @#if price_targeting ==1
                Pi=1;
            @#else
                @#if wage_targeting ==1
                    Pi_w=1;
                @#else
                    @#if composite_targeting ==1
                        Pi^(lambda_w/(lambda_w+lambda_p))*Pi_w^(lambda_p/(lambda_w+lambda_p))=1;
                    @#endif
                @#endif
            @#endif                
        @#else
            @#if price_targeting ==1
                R=1/betta*Pi^phi_pi;
            @#else
                @#if wage_targeting ==1
                    R=1/betta*Pi_w^phi_pi;
                @#else
                    @#if composite_targeting ==1
                        R=1/betta*(Pi^(lambda_w/(lambda_w+lambda_p))*Pi_w^(lambda_p/(lambda_w+lambda_p)))^phi_pi;
                    @#endif
                @#endif       
            @#endif
        @#endif
    @#endif
@#if Calvo
    [name='Market Clearing, eq. (15)']
    C=Y;
@#else
    [name='Market Clearing, eq. (15)']
    C=Y*(1-phi_w/2*(Pi_w/steady_state(Pi)-1)^2);
@#endif
    [name='Technology Shock, eq. (6)']
    log(A)=rho_a*log(A(-1))+eps_a;
    [name='Preference Shock, p.54']
    log(Z)=rho_z*log(Z(-1))-eps_z;
    [name='Monetary policy shock']
    nu=rho_nu*nu(-1)+eps_nu;
    [name='Definition marginal cost']
    MC=W_real*tau_s_p/((1-alppha)*Y/N_d*S);
    @#if sticky_prices == 1
        [name='LOM prices, eq. (7)']
        1=theta_p*Pi^(epsilon_p-1)+(1-theta_p)*(Pi_star)^(1-epsilon_p);
        [name='LOM price dispersion']
        S=(1-theta_p)*Pi_star^(-epsilon_p/(1-alppha))+theta_p*Pi^(epsilon_p/(1-alppha))*S(-1);
        [name='FOC price setting']
        Pi_star^(1+epsilon_p*(alppha/(1-alppha)))=x_aux_1/x_aux_2*epsilon_p/(epsilon_p-1);
        [name='Auxiliary price setting recursion 1']
        x_aux_1=Z*C^(-siggma)*Y*MC+betta*theta_p*Pi(+1)^(epsilon_p+alppha*epsilon_p/(1-alppha))*x_aux_1(+1);
        [name='Auxiliary price setting recursion 2']
        x_aux_2=Z*C^(-siggma)*Y+betta*theta_p*Pi(+1)^(epsilon_p-1)*x_aux_2(+1);
    @#else
        W_real*epsilon_p/(epsilon_p-1)*tau_s_p=(1-alppha)*A*N^(-alppha);
        S=1;
    @#endif
    
    [name='Definition log output']
    log_y = log(Y);
    [name='Definition log real wage']
    log_W_real=log(W_real);
    [name='Definition log hours']
    log_N=log(N);
    [name='Annualized inflation']
    pi_p_ann=4*log(Pi);
    [name='Annualized wage inflation']
    pi_w_ann=4*log(Pi_w);
    [name='Annualized nominal interest rate']
    i_ann=4*log(R);
    [name='Annualized real interest rate']
    r_real_ann=4*log(realinterest);
    [name='Real money demand, eq. (4)']
    M_real=Y/R^eta;
    [name='Definition log TFP']
    log_A=log(A);
    [name='Definition log preference']
    log_Pi=log(Pi);
    [name='Definition log preference']
    log_Pi_w=log(Pi_w);
    [name='Definition log preference']
    log_Z=log(Z);
    [name='Labor tax process']
    tau_n=(1-rho_tau_n)*tau_n_SS+rho_tau_n*tau_n(-1)+eps_tau_n;
    @#if Calvo
        @#if SGU_framework
            @#if log_utility
                [name='Period utility']
                Utility=Z*(log(C)-N^(1+varphi)/(1+varphi));
            @#else
                [name='Period utility']
                Utility=Z*(C^(1-siggma)/(1-siggma)-N^(1+varphi)/(1+varphi));
            @#endif
            [name='Auxiliary variable wage dispersion utility']        
            X_w=1;
        @#else
            @#if log_utility
                [name='Period utility']
                Utility=Z*(log(C)-N_d^(1+varphi)*X_w/(1+varphi));
            @#else
                [name='Period utility']
                Utility=Z*(C^(1-siggma)/(1-siggma)-N_d^(1+varphi)*X_w/(1+varphi));
            @#endif
            [name='Auxiliary variable wage dispersion utility']        
            X_w=(1-theta_w)*Pi_w_star^(-epsilon_w*(1+varphi))+theta_w*Pi_w^(epsilon_w*(1+varphi))*X_w(-1);
        @#endif
    @#else
        @#if log_utility
            [name='Period utility']
            Utility=Z*(log(C)-N^(1+varphi)/(1+varphi));
        @#else
            [name='Period utility']
            Utility=Z*(C^(1-siggma)/(1-siggma)-N^(1+varphi)/(1+varphi));
        @#endif
        [name='Auxiliary variable wage dispersion utility']        
        X_w=1;
    @#endif
    @#if Ramsey_policy==0    
        [name='Recursive Welfare']
        Recursive_Welfare=Utility+betta*Recursive_Welfare(+1);
    @#endif
    @#if log_utility
        [name='Natural Recursive Welfare']
        Recursive_natural_welfare_equivalent=Z*(log((1-lambda_utility)*Y_nat)-N_nat^(1+varphi)/(1+varphi))+betta*Recursive_natural_welfare_equivalent(+1);
    @#else
        [name='Natural Recursive Welfare']
        Recursive_natural_welfare_equivalent=Z*(((1-lambda_utility)*Y_nat)^(1-siggma)/(1-siggma)-N_nat^(1+varphi)/(1+varphi))+betta*Recursive_natural_welfare_equivalent(+1);
    @#endif
    @#if Ramsey_policy==0    
        [name='Definition welfare gap']
        Welfare_gap=Recursive_Welfare-Recursive_natural_welfare_equivalent;
    @#endif        
    [name='Definition marginal product of labor']
    MPN=(1-alppha)*A*(N_d/S)^(1-alppha)/N_d;
    [name='Definition after tax marginal rate of substitution (not sure about correct aggregation/concept of aggregate MRS, only used to check steady state)']
    MRS=(1+tau_c)/(1+tau_n)*N^varphi/C^(-siggma);
    [name='Definition natural hours worked (MRS=MPN with Y=C)']                
    (1+tau_c)/(1+tau_n)*N_nat^varphi/(A*N_nat^(1-alppha))^(-siggma)*epsilon_w/(epsilon_w-1)/tau_s_w
            =(1-alppha)*A*N_nat^(-alppha)*(epsilon_p-1)/epsilon_p/tau_s_p;
    [name='Definition natural output']            
    Y_nat=A*N_nat^(1-alppha);
    [name='Definition log output gap']
    log_Y_gap=log(Y)-log(Y_nat);    
end;

%----------------------------------------------------------------
%  Steady state values
%---------------------------------------------------------------

steady_state_model;
tau_n=tau_n_SS;
lambda_p=(1-theta_p)*(1-betta*theta_p)/theta_p*((1-alppha)/(1-alppha+alppha*epsilon_p));      %defined on page 166, Gali (2015)
@#if fixed_WPC_slope
    %keep fixed slope of WPC and compute parameters accordingly
    @#if SGU_framework
        @#if Calvo
            theta_w=get_Calvo_theta(lambda_w,epsilon_w,betta,varphi,1);
        @#else
            phi_w=(epsilon_w-1)*(1-tau_n_SS)/lambda_w*(1-alppha)*(epsilon_p-1)/epsilon_p/tau_s_p*tau_s_w;     
        @#endif
    @#else
        @#if Calvo
            theta_w=get_Calvo_theta(lambda_w,epsilon_w,betta,varphi,0);
        @#else
            phi_w=(epsilon_w-1)*(1-tau_n_SS)/lambda_w*(1-alppha)*(epsilon_p-1)/epsilon_p/tau_s_p*tau_s_w;     
        @#endif
    @#endif
@#else
    @#if SGU_framework
        %compute slope (and ) based on Calvo parameter 
        @#if Calvo
            lambda_w=(1-theta_w)*(1-betta*theta_w)/theta_w;
        @#else
            phi_w=(epsilon_w-1)/((1-theta_w)*(1-betta*theta_w))*(1-tau_n)*theta_w*(1-alppha)*(epsilon_p-1)/epsilon_p/tau_s_p*tau_s_w;
            lambda_w=(epsilon_w-1)/phi_w*(1-tau_n)*(1-alppha)*(epsilon_p-1)/epsilon_p/tau_s_p*tau_s_w;
        @#endif
    @#else
        @#if Calvo
            lambda_w=(1-theta_w)*(1-betta*theta_w)/(theta_w*(1+epsilon_w*varphi));
        @#else
            phi_w=(epsilon_w-1)/((1-theta_w)*(1-betta*theta_w))*(1-tau_n)*theta_w*(1-alppha)*(epsilon_p-1)/epsilon_p/tau_s_p*tau_s_w*(1+epsilon_w*varphi);
            lambda_w=(epsilon_w-1)*(1-tau_n)/phi_w*(1-alppha)*(epsilon_p-1)/epsilon_p/tau_s_p*tau_s_w;
        @#endif
    @#endif
@#endif
A=1;
Z=1;
S=1;
S_w=1;
Pi_star=1;
@#if Calvo
    Pi_w_star=1;        
@#endif
MC=(epsilon_p-1)/epsilon_p;
@#if Ramsey_policy==0
    R=1/betta;
@#endif
Q=1/R;
realinterest=1/betta;
Pi=R/realinterest;
Pi_w=Pi;            
N=((epsilon_w-1)/epsilon_w*(tau_s_w/tau_s_p)*(1-tau_n)/(1+tau_c)*(1-alppha)*MC)^(1/((1-siggma)*alppha+varphi+siggma));
C=A*N^(1-alppha);
W_real=epsilon_w/(epsilon_w-1)/tau_s_w*C^siggma*N^varphi*(1+tau_c)/(1-tau_n);
Y=C;
nu=0;
x_aux_1=C^(-siggma)*Y*MC/(1-betta*theta_p*Pi^(epsilon_p/(1-alppha)));
x_aux_2=C^(-siggma)*Y/(1-betta*theta_p*Pi^(epsilon_p-1));
N_d=N/S_w;    
@#if Calvo
    @#if SGU_framework
        f=(epsilon_w-1)/epsilon_w*tau_s_w*(1-tau_n)/(1+tau_c)*C^(-siggma)*W_real*Pi_w_star^(1-epsilon_w)*N_d/
                    (1-betta*theta_w*(Pi_w)^(epsilon_w-1));
    @#else
        f=(epsilon_w-1)/epsilon_w*tau_s_w*(1-tau_n)/(1+tau_c)*C^(-siggma)*W_real*Pi_w_star^(-epsilon_w)*N_d/
                    (1-betta*theta_w*(Pi_w)^(epsilon_w-1));
    @#endif
@#endif
log_y = log(Y);
log_W_real=log(W_real);
log_N=log(N);
pi_p_ann=4*log(Pi);
pi_w_ann=4*log(Pi_w);
i_ann=4*log(R);
r_real_ann=4*log(realinterest);
M_real=Y/R^eta;
log_A=0;
log_Z=0;
log_Pi=log(Pi);
log_Pi_w=log(Pi_w);
MPN=(1-alppha)*A*(N_d/S)^(1-alppha)/N_d;
MRS=(1+tau_c)/(1+tau_n)*N^varphi/C^(-siggma);
N_nat=N;
Y_nat=Y;
@#if log_utility
    Utility=(log(C)-N^(1+varphi)/(1+varphi));
    @#if Ramsey_policy==0    
        Recursive_Welfare=1/(1-betta)*(log(C)-N^(1+varphi)/(1+varphi));
    @#endif    
    Recursive_natural_welfare_equivalent=1/(1-betta)*(log((1-lambda_utility)*Y_nat)-N_nat^(1+varphi)/(1+varphi));
@#else
    Utility=(C^(1-siggma)/(1-siggma)-N^(1+varphi)/(1+varphi));
    @#if Ramsey_policy==0    
        Recursive_Welfare=1/(1-betta)*Utility;
    @#endif    
    Recursive_natural_welfare_equivalent=1/(1-betta)*(((1-lambda_utility)*Y_nat)^(1-siggma)/(1-siggma)-N_nat^(1+varphi)/(1+varphi));
@#endif
@#if Ramsey_policy==0    
    Welfare_gap=Recursive_Welfare-Recursive_natural_welfare_equivalent;
@#endif    
X_w=1;
log_Y_gap=0;    
end;

write_latex_original_model;//(write_equation_tags);
write_latex_static_model;
//write_latex_steady_state_model;
// write_latex_parameter_table;

% resid(1);
% steady;
% check;

%----------------------------------------------------------------
%  Do demand shocks
%---------------------------------------------------------------

shocks;
    var eps_a       = 0.01^2; 
    var eps_z       = 0; 
end;

@#if Ramsey_policy==1
    planner_objective Utility;
    ramsey_model(instruments=(R),planner_discount=betta);
    initval;
    R=1/betta;
end;
@#endif

stoch_simul(order=2,irf=0) log_Pi log_Pi_w log_Y_gap log_N;

%read out variable positions
pi_p_pos=strmatch('log_Pi',var_list_ ,'exact');
pi_w_pos=strmatch('log_Pi_w',var_list_ ,'exact');
y_gap_pos=strmatch('log_Y_gap',var_list_ ,'exact');

%read out variances
variance.pi_p=oo_.var(pi_p_pos,pi_p_pos);
variance.pi_w=oo_.var(pi_w_pos,pi_w_pos);
variance.y_gap=oo_.var(y_gap_pos,y_gap_pos);

%compute consumption equivalent
options_old=options_;
options_.nocorr=1;
options_.noprint=1;
lambda_unconditional_technology=csolve('get_consumption_equivalent_unconditional_welfare',0,[],1e-8,1000)
lambda_conditional_technology=csolve('get_consumption_equivalent_conditional_welfare',lambda_unconditional_technology,[],1e-8,1000)
options_=options_old;
 
%display results
labels={'sigma(pi_p)';'sigma(pi_w)';'sigma(tilde y)';'L unc';'L cond'};
headers={'Technology shocks';case_title};
values_technology=[sqrt([variance.pi_p;variance.pi_w;variance.y_gap]);lambda_unconditional_technology;lambda_conditional_technology];
options_.noprint=0;
dyntable(options_,table_title,headers,labels,100*values_technology,size(labels,2)+2,4,3)
%----------------------------------------------------------------
%  Do demand shocks
%---------------------------------------------------------------

shocks;
    var eps_z       = 0.01^2; 
    var eps_a       = 0; 
end;

stoch_simul(order=2,irf=0) log_Pi log_Pi_w log_Y_gap log_N log_y;

%read out variable positions
pi_p_pos=strmatch('log_Pi',var_list_ ,'exact');
pi_w_pos=strmatch('log_Pi_w',var_list_ ,'exact');
y_gap_pos=strmatch('log_Y_gap',var_list_ ,'exact');

%read out variances
variance.pi_p=oo_.var(pi_p_pos,pi_p_pos);
variance.pi_w=oo_.var(pi_w_pos,pi_w_pos);
variance.y_gap=oo_.var(y_gap_pos,y_gap_pos);

%compute consumption equivalent
options_old=options_;
options_.nocorr=1;
options_.noprint=1;
lambda_unconditional_demand=csolve('get_consumption_equivalent_unconditional_welfare',0,[],1e-8,1000)
lambda_conditional_demand=csolve('get_consumption_equivalent_conditional_welfare',lambda_unconditional_demand,[],1e-8,1000)
options_=options_old;

%display results
labels={'sigma(pi_p)';'sigma(pi_w)';'sigma(tilde y)';'L unc';'L cond'};
headers={'Demand shocks   ';case_title};
values_demand=[sqrt([variance.pi_p;variance.pi_w;variance.y_gap]);lambda_unconditional_demand;lambda_conditional_demand];
options_.noprint=0;
dyntable(options_,table_title,headers,labels,100*values_demand,size(labels,2)+2,5,4)      