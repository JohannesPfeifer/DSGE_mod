/*
 * This file replicates the Baseline model under Financial Autarky of
 * Ghironi/Melitz (2005), "International trade and macroeconomic dynamics with 
 * heterogeneous firms", Quarterly Journal of Economics, 120(3), 865-915
 * 
 * Notes: 
 *  - The file replicates Figure III on page 892 
 *  - Fabio Ghironi has pointed out three issues with the original paper at http://faculty.washington.edu/ghiro/GhiroMelitzQJE05Bugs.pdf
 *      Point 1 affects the present file. The formula used to set the parameter k is wrong. Thus, the actual calibration 
 *      target is not satisfied.
 *  - The loglinear option is used to get percentage deviations from steady state, 
 *    while the model is entered in levels
 *  - The steady_state-file is used to handle parameter dependence, i.e. to set 
 *    the fixed costs as a share of flow entry costs
 *  - The interest r_t in the paper is a predetermined variable as it is known at t-1; 
 *    the timing has been adjusted to reflect Dynare's timing convention
 *  - A detailed derivation of the condition for Average productivities can be found in
 *    the Derivation_average_productivities.pdf on Github
 *  - The Latex names use a workaround for tildes. See https://tex.stackexchange.com/questions/32739/double-superscript-error-involving-tilde
 *
 * This file was originally written by William Gatt, University of Nottingham & Central Bank of Malta.
 * It has been updated by Johannes Pfeifer and Lucas Radke. In case you spot mistakes, 
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model
 */

/*
 * Copyright (C) 2017 William Gatt
 *               2018 Johannes Pfeifer and Lucas Radke
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Endogeneous variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var ztilded         ${{}\tilde z_{D}}$                (long_name='Average productivity level - home economy')
    ztilded_        ${{}\tilde z^{\star}_{D}}$        (long_name='Average productivity level - foreign economy')
    ztildex         ${{}\tilde z_{X}}$                (long_name='Average productivity level exporters - home economy') 
    ztildex_        ${{}\tilde z^{\star}_{X}}$        (long_name='Average productivity level exporters - foreign economy')
    zx              ${z_{X}}$                         (long_name='Productivity cutoff level - home economy')    
    zx_             ${z^{\star}_{X}}$                 (long_name='Productivity cutoff level - foreign economy')
    w               ${w}$                             (long_name='Real wage - home economy') 
    w_              ${w^{\star}}$                     (long_name='Real wage - foreign economy')
    rhotilded       ${{}\tilde \rho_{D}}$             (long_name='Average relative prices of home producers in the home market') 
    rhotilded_      ${{}\tilde \rho^{\star}_{D}}$     (long_name='Average relative prices of foreign producers in the foreign market')
    rhotildex       ${{}\tilde \rho_{X}}$             (long_name='Average relative prices of foreign exporters in the home market')
    rhotildex_      ${{}\tilde \rho^{\star}_{X}}$     (long_name='Average relative prices of home exporters in the foreign market')
    dtilde          ${{}\tilde d}$                    (long_name='Average total profits of home and foreign firms - home economy')
    dtilde_         ${{}\tilde d^{\star}}$            (long_name='Average total profits of home and foreign firms - foreign economy')
    dtilded         ${{}\tilde d_{D}}$                (long_name='Average firm profit earned from domestic sales for all home producers')
    dtilded_        ${{}\tilde d^{\star}_{D}}$        (long_name='Average firm profit earned from domestic sales for all foreign producers')
    dtildex         ${{}\tilde d_{X}}$                (long_name='Average firm export profits for all home exporters')
    dtildex_        ${{}\tilde d^{\star}_{X}}$        (long_name='Average firm export profits for all foreign exporters')
    Ne              ${N_{E}}$                         (long_name='Mass of entrants - home economy')
    Ne_             ${N^{\star}_{E}}$                 (long_name='Mass of entrants - foreign economy')
    Nd              ${N_{D}}$                         (long_name='Number of home-producing firms - home economy')
    Nd_             ${N^{\star}_{D}}$                 (long_name='Number of home-producing firms - foreign economy')
    Nx              ${N_{X}}$                         (long_name='Number of home firms exporting to foreign market')
    Nx_             ${N^{\star}_{X}}$                 (long_name='Number of foreign firms exporting to home market')
    r               ${r}$                             (long_name='Interest rate on holdings of bonds - home economy')
    r_              ${r^{\star}}$                     (long_name='Interest rate on holdings of bonds - home economy')
    vtilde          ${{}\tilde v}$                    (long_name='Present discounted value of home firms expected stream of profits')  
    vtilde_         ${{}\tilde v^{\star}}$            (long_name='Present discounted value of home firms expected stream of profits') 
    C               ${C}$                             (long_name='Consumption - home economy')
    C_              ${C^{\star}}$                     (long_name='Consumption - foreign economy')
    Z               ${Z}$                             (long_name='Aggregate labor productivity - home economy')
    Z_              ${Z^{\star}}$                     (long_name='Aggregate labor productivity - foreign economy')
    Q               ${Q}$                             (long_name='Real exchange rate - welfare based')
    Qtilde          ${{}\tilde Q}$                    (long_name='Real exchange rate - empirical')
    TOL             ${TOL}$                           (long_name='Terms of labor')
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exogenous  variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varexo eps_Z           ${\varepsilon_{Z}}$               (long_name='Productivity shock - home economy')
       eps_Z_          ${\varepsilon^{\star}_{Z}}$       (long_name='Productivity shock - foreign economy')                        
;


parameters
    theta       ${\theta}$             (long_name='elasticity of subsitution across goods')
    k           ${k}$                  (long_name='shape parameter')
    delta       ${\delta}$             (long_name='exogenous exit rate')
    gamma       ${\gamma}$             (long_name='relative risk aversion')
    betaa       ${\beta}$              (long_name='discount factor')
    tau         ${\tau}$               (long_name='iceberg cost')
    zmin        ${z_{min}}$            (long_name='minimum relative prod - home')
    zmin_       ${z^{\star}_{min}}$    (long_name='minimum relative prod - foreign')
    fe          ${f_{E}}$              (long_name='entry cost - home')
    fe_         ${f^{\star}_{E}}$      (long_name='entry cost - foreign')
    fx          ${f_{X}}$              (long_name='fixed cost of exporting - home')
    fx_         ${f^{\star}_{X}}$      (long_name='fixed cost of exporting - foreign')
    L           ${L}$                  (long_name='Labour supply - home')
    L_          ${L^{\star}}$          (long_name='Labour supply - foreign')
    rhoZ        ${\rho_{Z}}$            (long_name='Productivity AR process parameter - home')
    rhoZ_       ${\rho^{\star}_{Z}}$    (long_name='Productivity AR process parameter - foreign')
    fx_share    ${FCshare}$             (long_name='Share of fixed cost as share of flow value of entry cost')
   ;
    
    % Calibration
    betaa   =   0.99;
    gamma   =   2.0 ;
    delta   =   0.025;
    theta   =   3.8 ; 
    k       =   3.4 ; %see note above
    tau     =   1.3 ;
    zmin    =   1.0 ;
    zmin_   =   1.0 ;
    fe      =   1.0 ; 
    fe_     =   1.0 ; 
    L       =   1.0 ; 
    L_      =   1.0 ;
    rhoZ    =   0.9 ;
    rhoZ_   =   0.9 ;
    fx_share=   0.235;

model;
%[name='constant linking average productivites']
# vv = (k/(k-(theta-1)))^(1/(theta-1));     
% Price indices
  [name='Price indices - home economy']
   1            =  Nd*(rhotilded)^(1-theta) + Nx_*(rhotildex_)^(1-theta);
  [name='Price indices - foreign economy']
   1            =  Nd_*(rhotilded_)^(1-theta) + Nx*(rhotildex)^(1-theta);
  [name='Average real domestic price - home economy']
   rhotilded    = (theta/(theta-1))*(w/(Z*ztilded));
  [name='Average real domestic price - foreign economy']
   rhotilded_   = (theta/(theta-1))*(w_/(Z_*ztilded_));
  [name='Average real export price - home economy'] 
   rhotildex    = (tau*(theta/(theta-1))*(w/(Z*ztildex)))/Q;
  [name='Average real export price - foreign economy']
   rhotildex_   = (tau*(theta/(theta-1))*(w_/(Z_*ztildex_)))*Q;
   
% Profits
   [name='Average total profits of home firms']
    dtilde      = dtilded  + (Nx/Nd)*dtildex;
   [name='Average total profits of foreign firms']
    dtilde_     = dtilded_ + (Nx_/Nd_)*dtildex_;
   [name='Average profits earned from domestic sales - home economy'] 
    dtilded     = (1/theta)*rhotilded^(1-theta)*C;
   [name='Average profits earned from domestic sales - foreign economy']
    dtilded_    = (1/theta)*rhotilded_^(1-theta)*C_;
    
% Free entry
   [name='Free entry condition - home economy']
    vtilde      = w*(fe/Z);
   [name='Free entry condition - foreign economy']
    vtilde_     = w_*(fe_/Z_);
    
% Zero-profit export cutoffs
   [name='Zero-profit export cutoffs - home economy']
    dtildex     =  w*(fx/Z)*((theta-1)/(k-(theta-1)));
   [name='Zero-profit export cutoffs - foreign economy']
    dtildex_    =  w_*(fx_/Z_)*((theta-1)/(k-(theta-1)));

% Share of exporting firms
   [name='Share of exporting firms - home economy']
    Nx/Nd       =  ((zmin/ztildex)^k)*((k/(k-(theta-1)))^(k/(theta-1)));
   [name='Share of exporting firms - foreign economy']
    Nx_/Nd_     =  ((zmin_/ztildex_)^k)*((k/(k-(theta-1)))^(k/(theta-1)));

% Number of firms
   [name='Number of firms - home economy']
    Nd          = (1-delta)*(Nd(-1) + Ne(-1));
   [name='Number of firms - home economy']
    Nd_         = (1-delta)*(Nd_(-1) + Ne_(-1));

% Euler equation (bonds)
   [name='Euler equation (bonds) - home economy']
    C^(-gamma)  =   betaa*(1+r)*(C(+1)^(-gamma));
   [name='Euler equation (bonds) - foreign economy']
    C_^(-gamma) =   betaa*(1+r_)*(C_(+1)^(-gamma));
    
% Euler equation (shares)
   [name='Euler equation (shares) - home economy']
    vtilde      =   betaa*(1-delta)*(((C(+1)/C)^(-gamma))*(vtilde(+1)+ dtilde(+1)));
   [name='Euler equation (shares) - foreign economy']
    vtilde_     =   betaa*(1-delta)*(((C_(+1)/C_)^(-gamma))*(vtilde_(+1)+ dtilde_(+1)));

% Aggregate accounting
   [name='Aggregate accounting - home economy']
    C           = w*L + Nd*dtilde - Ne*vtilde;
   [name='Aggregate accounting - foreign economy']
    C_          = w_*L_ + Nd_*dtilde_ - Ne_*vtilde_;
    
% Balanced trade (real exchange rate - welfare based)
   [name='Real exchange rate - welfare based']
    Q           = (Nx_*(rhotildex_^(1-theta))*C)/(Nx*(rhotildex^(1-theta))*C_);
    
% Real exchange rate - empirical 
   [name='Definition terms of labor (eq. (4))']
    Qtilde      = (((Nd_/(Nd_+Nx))*((TOL^(1-theta))) + (Nx/(Nd_+Nx))*((tau*(ztilded/ztildex))^(1-theta)))/((Nd/(Nd+Nx_))+((Nx_/(Nd+Nx_))*((TOL*tau*(ztilded_/ztildex_))^(1-theta)))))^(1/(1-theta));
   [name='Link between empirical and theoretical exchange rate measures (p. 880)']
    Qtilde      = Q*((Nd + Nx_)/(Nd_+Nx))^(-1/(theta-1));
    
% Aggregate productivity
   [name='Aggregate productivity - home economy']
    Z           =   (1-rhoZ)*1.0 + rhoZ*Z(-1) + eps_Z;
   [name='Aggregate productivity - foreign economy']
    Z_          =   (1-rhoZ_)*1.0 + rhoZ_*Z_(-1) + eps_Z_;
    
% Average producivity
   [name='Relation average and cutoff productivity - home economy']
    ztilded     = vv*zmin;
   [name='Relation average and cutoff productivity - foreign economy']
    ztilded_    = vv*zmin_;
    
   [name='Average productivity - home economy']
    ztildex     = ((w/Z)^(theta)*fx*((theta-1)/(k-(theta-1))+1)*(Q^(-theta))*(tau^(theta-1))*theta*((theta/(theta-1))^(theta-1))*(C_^(-1)))^(1/(theta-1));
   [name='Average productivity - foreign economy']
    ztildex_    = ((w_/Z_)^(theta)*fx_*((theta-1)/(k-(theta-1))+1)*(Q^(theta))*(tau^(theta-1))*theta*((theta/(theta-1))^(theta-1))*(C^(-1)))^(1/(theta-1));

   [name='Export cut-off productivities - home economy']
    zx          = ztildex/vv;
   [name='Export cut-off productivities - foreign economy']
    zx_         = ztildex_/vv;

end;

// options_.TeX=1;

% Shock to productivity
shocks;
    var eps_Z = 0.0001;
end;

steady;
check;

stoch_simul(order=1,nograph,hp_filter=1600,irf=100,loglinear) C C_ Ne Nd zx Nx Ne_ Nd_ zx_ Nx_ TOL Qtilde Q Z ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATING FIGURE III ON PAGE 892 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
figure('Name','Replication of Figure III');
Horizon = (1:length(oo_.irfs.C_eps_Z)+1)/4;
zeroline = (zeros(length(oo_.irfs.C_eps_Z)+1,1))';

subplot(4,4,1)
plot(Horizon,[0 , (oo_.irfs.C_eps_Z)*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('C','Fontsize',12);
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([0 0.5])

subplot(4,4,2)
plot(Horizon,[0 , (oo_.irfs.C__eps_Z)*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('C*','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.05 0.05]);

subplot(4,4,3)
plot(Horizon,[0 , (oo_.irfs.Ne_eps_Z)*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('Ne','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.25 4]);

subplot(4,4,4)
plot(Horizon,[0 , oo_.irfs.Nd_eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('Nd','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.05 0.55]); 

subplot(4,4,5)
plot(Horizon,[0 , oo_.irfs.zx_eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('z_{x}','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.05 0.17]);

subplot(4,4,6)
plot(Horizon,[0 , oo_.irfs.Nx_eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('Nx','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.1 0.2]);

subplot(4,4,7)
plot(Horizon,[0 , oo_.irfs.Ne__eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('Ne*','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.05 0.05]);

subplot(4,4,8)
plot(Horizon,[0 , oo_.irfs.Nd__eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('Nd*','Fontsize',12)    
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.05 0.05]); 

subplot(4,4,9)
plot(Horizon,[0 , oo_.irfs.zx__eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('z_{x}*','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.05 0.05]);

subplot(4,4,10)
plot(Horizon,[0 , oo_.irfs.Nx__eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('Nx*','Fontsize',12)        
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([0 0.15]);

subplot(4,4,11)
plot(Horizon,[0 , oo_.irfs.TOL_eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('TOL','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.12 0.07]);

subplot(4,4,12)
plot(Horizon,[0 , oo_.irfs.Qtilde_eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('Q~','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.12 0.05]);

subplot(4,4,13)
plot(Horizon,[0 , oo_.irfs.Q_eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('Q','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([-0.05 0.055]);

subplot(4,4,14)
plot(Horizon,[0 , oo_.irfs.Z_eps_Z*100],'Color','r','LineWidth',2.5,'LineStyle','-'); hold on
plot(Horizon,zeroline,'Color','k','LineWidth',1.0);
title('Z','Fontsize',12)
ylabel('% deviation from SS', 'FontSize', 8)
xlim([0 max(Horizon)])
ylim([0 1]);
   
// %Latex-Configuration
// write_latex_original_model;
// write_latex_static_model;
// write_latex_dynamic_model(write_equation_tags);
// write_latex_parameter_table;
// write_latex_definitions;
// collect_latex_files;