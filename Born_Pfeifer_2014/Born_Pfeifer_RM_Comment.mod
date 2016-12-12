/*
 * This file allows replicating the results of:
 * Benjamin Born and Johannes Pfeifer (2014): "Risk Matters: A comment", American Economic Review
 * 
 * It provides a full replication of the results of the paper.
 * 
 * Notes:
 * - This mod-file requires at least Dynare 4.4.0 and has been tested with Dynare 4.4.3 using Matlab 2013a and 2014b
 * - It uses the Matlab Econometrics Toolbox for HP-filtering
 * - For comparability, we follow Jesús Fernández-Villaverde, Pablo Guerrón-Quintana, Juan F. Rubio-Ramírez, 
 *      and Martin Uribe (2011): "Risk Matters" in their simulation approach, American Economic Review 101 
 *       (October 2011): 2530-2561.
 *   This in particular implies that:
 *      - Their pruning is used (implemented by Johannes Pfeifer in simult_FGRU.m)
 *      - The first period shock vector is always the same across all replications
 *      - The first entry of the simulated series is always the deterministic steady state
 *      - The last simulated entry for the first order shock term functions as the starting value of the next simulation run
 * - This pruning scheme is non-standard. To use the standard Andreasen et al. (2013)-pruning scheme, replace the calls to
 *      simult_FGRU.m with calls to simult_.m from Dynare (see the comments below)
 * - The runtime with 10000 simulation replications can be between 15 and 30 Minutes, depending on your computer
 *
 *
 * Copyright (C) 2013-14 Benjamin Born and Johannes Pfeifer
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
 * For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.
 */

//define number of simulation repetitions
replications=10000;
estimation_replications=200;

//set country; you must not set more than one entry to 1
@#define Argentina = 1
@#define Ecuador = 0
@#define Venezuela = 0
@#define Brazil = 0
    
//set this to 1 for the recalibrated model and to 0 to use the original calibration      
@#define recalibration = 0
        
//set this to 1 to estimate the model after running the simulations
@#define estimation = 0



//define strings for save name
@#if recalibration == 1
    end_save_string='recalibration';
@# else
    end_save_string='replication';
@# endif

//******************************Dynare definitions************************************************

var sigma_r sigma_tb eps_r eps_tb X  D K lambda C H Y  I phi r NX NX_Y CA CA_Y NX_Y_quarterly;
varexo u_x u_tb u_r u_sigma_tb u_sigma_r;

predetermined_variables K D; //define predetermined variables

parameters r_bar rho_eps_r sigma_r_bar rho_sigma_r eta_r 
           rho_eps_tb sigma_tb_bar rho_sigma_tb eta_tb
            delta alppha nu rho_x betta
            Phi phipar sigma_x D_bar thetheta eta;

//******************************Dynare parameter initialization************************************************

// Calibration from Table 4 of FGRU
rho_eps_tb=0.95;
sigma_tb_bar=-8.06;     //8.05 in paper, but 8.06 in code
rho_sigma_tb=0.94;
eta_tb=0.13;
nu=5;                   //inverse of the elasticity of intertemporal substitution
eta=1000;               //elasticity of labor to wages
delta=0.014;            //depreciation
alppha = 0.32;          //capital income share
rho_x=0.95;             //autocorrelation TFP


@#if Argentina == 1 
    country_string='Argentina';
    // Calibration from Table 3
    rho_eps_r=0.97;
    sigma_r_bar=-5.71;
    rho_sigma_r=0.94; 
    eta_r=0.46;

    r_bar=log(0.02);
    betta = 1/(1+exp(r_bar));   //discount factor

    @#if recalibration == 0 
        //original calibration
        Phi=0.001; %debt elasticity
        D_bar= 4; % steady state debt
        phipar=95; % capital adjustment costs
        sigma_x=log(0.015); 
    @# else
        //Recalibration from Born/Pfeifer's Table 
        sigma_x=-3.21695494786718;
        phipar=47.8376372268423;
        D_bar= 18.8015551328765;
        Phi=log(1.00059258914564);
    @# endif
@# endif

@#if Ecuador == 1 
    country_string='Ecuador';
    // Calibration from Table 3
    rho_eps_r=0.95;
    sigma_r_bar=-6.06;
    rho_sigma_r=0.96; 
    eta_r=0.35;

    // Calibration from Table 6 //
    r_bar=log(0.011);
    betta = 1/(1+exp(r_bar)); %discount factor

    @#if recalibration == 0 
        //original calibration
        Phi=0.001; %debt elasticity
        D_bar= 13; % steady state debt
        phipar=35; % capital adjustment costs
        sigma_x=log(0.0055); 
    @# else
        sigma_x=-4.48184679235884; 
        phipar=18.9815814953814; % capital adjustment costs
        D_bar= 26.3037546297346; % steady state debt
        Phi=log(1.00015105501358); %debt elasticity
    @# endif
@# endif


@#if Venezuela == 1 
    country_string='Venezuela';
    // Calibration from Table 3
    rho_eps_r=0.94;
    sigma_r_bar=-6.88;
    rho_sigma_r=0.91; 
    eta_r=0.32;

    // Calibration from Table 6 //
    r_bar=log(0.007);
    betta = 1/(1+exp(r_bar)); %discount factor

    @#if recalibration == 0 
        //original calibration
        Phi=5e-4; %debt elasticity
        D_bar= 22; % steady state debt
        phipar=12; % capital adjustment costs
        sigma_x=log(0.013); 
    @# else
        sigma_x=-3.19414337673133; 
        phipar=5.08058763451490; % capital adjustment costs
        D_bar= 21.4624425882826; % steady state debt
        Phi=log(22026.4657954966); %debt elasticity
    @# endif
@# endif


@#if Brazil == 1 
    country_string='Brazil';
    // Calibration from Table 3
    rho_eps_r=0.95;
    sigma_r_bar=-6.97;
    rho_sigma_r=0.95; 
    eta_r=0.28;

    // Calibration from Table 6 //
    r_bar=log(0.007);
    betta = 1/(1+exp(r_bar)); %discount factor

    @#if recalibration == 0 
        //original calibration
        Phi=5e-4; %debt elasticity
        D_bar= 3; % steady state debt
        phipar=50; % capital adjustment costs
        sigma_x=log(0.013); 
    @# else
        sigma_x=-3.18192192158336; 
        phipar=101.055881220935; % capital adjustment costs
        D_bar= 0.553728179249502; % steady state debt
        Phi=log(5.40237854463379); %debt elasticity
    @# endif
@# endif



thetheta=1; //normalization parameter from FGRU; unspecified in their manuscript and extracted from their code

//******************************Dynare model equations************************************************
model;
(exp(C))^-nu=exp(lambda);
exp(lambda)/(1+r)=exp(lambda)*Phi*((D(+1))-D_bar)+betta*exp(lambda(+1));
-exp(phi)+betta*((1-delta)*exp(phi(+1))+alppha*exp(Y(+1))/exp(K(+1))*exp(lambda(+1)))=0;
thetheta*exp(H)^eta=(1-alppha)*exp(Y)/exp(H)*exp(lambda);
exp(phi)*(1-phipar/2*((exp(I)-exp(I(-1)))/exp(I(-1)))^2-phipar*exp(I)/exp(I(-1))*((exp(I)-exp(I(-1)))/exp(I(-1))))+betta*exp(phi(+1))*phipar*(exp(I(+1))/exp(I))^2*((exp(I(+1))-exp(I))/exp(I))=exp(lambda);

exp(Y)=exp(K)^alppha*(exp(X)*exp(H))^(1-alppha);
X=rho_x*X(-1)+exp(sigma_x)*u_x;
exp(K(+1))=(1-delta)*exp(K)+(1-phipar/2*(exp(I)/exp(I(-1))-1)^2)*exp(I);
NX=exp(Y)-exp(C)-exp(I);
NX_Y=NX/exp(Y);
CA=D-D(+1);
CA_Y=CA/exp(Y);
NX_Y_quarterly=(NX+NX(-1)+NX(-2))/(exp(Y)+exp(Y(-1))+exp(Y(-2)));
exp(Y)-exp(C)-exp(I)=(D)-(D(+1))/(1+r)+Phi/2*((D(+1))-D_bar)^2;
r=exp(r_bar)+eps_tb+eps_r;
eps_tb=rho_eps_tb*eps_tb(-1)+exp(sigma_tb)*u_tb;   
sigma_tb=(1-rho_sigma_tb)*sigma_tb_bar+rho_sigma_tb*sigma_tb(-1)+eta_tb*u_sigma_tb;
eps_r=rho_eps_r*eps_r(-1)+exp(sigma_r)*u_r;   
sigma_r=(1-rho_sigma_r)*sigma_r_bar+rho_sigma_r*sigma_r(-1)+eta_r*u_sigma_r;
end;


//******************************Initial values for steady state computation******************************************
initval;
sigma_tb=sigma_tb_bar;
sigma_r=sigma_r_bar;
eps_r=0;
eps_tb=0;
D=D_bar;

% steady states taken from FGRU Mathematica code
C=0.8779486025329908;
K=3.293280327636415;
lambda=-4.389743012664954;
H=-0.0037203652717462993;
phi=-4.389743012664954;
I=-0.9754176217303792;
Y=1.0513198564588924;
r=exp(r_bar);
NX=exp(Y)-exp(C)-exp(I); 
NX_Y=NX/exp(Y);
NX_Y_quarterly=NX_Y;
CA=0;
CA_Y=0;
end;

//******************************Set covariance matrix of structural shocks to identity**************************************
shocks;
var u_x; stderr 1;
var u_r; stderr 1;
var u_tb; stderr 1;
var u_sigma_tb; stderr 1;
var u_sigma_r; stderr 1;
end;

//******************************Compute steady state**************************************
options_.solve_tolf=1E-12;
steady(solve_algo=3);
check;

//get decision rules at order 3;
stoch_simul(order=3,pruning,periods=0,irf=0,nofunctions); 



npts=96;

//******************************get ergodic mean in the absence of shocks as baseline**************************************
//******************************Note: we don't do an additional convergence check here*************************************
//******************************Other parametrizations might require it+++++++++++++++*************************************

burnin=4000; //periods for convergence
shocks_mat_with_zeros=zeros(burnin,M_.exo_nbr); //shocks set to 0 to simulate without uncertainty
%% Note that simult_FGRU.m uses the original non-standard pruning scheme
out_noshock = simult_FGRU(oo_.dr.ys,oo_.dr,shocks_mat_with_zeros,options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1)); //simulate series
%% To simulate serires with proper pruning scheme as in Andreasen et al. (2013), use the following code
% out_noshock = simult_(oo_.dr.ys,oo_.dr,shocks_mat_with_zeros,options_.order); //simulate series

log_deviations_SS_noshock=out_noshock-oo_.dr.ys*ones(1,burnin+M_.maximum_lag); //subtract steady state to get deviations from steady state
ergodicmean_no_shocks=out_noshock(:,end); //EMAS is the final point
nx_y_EMAS=ergodicmean_no_shocks(strmatch('NX_Y',M_.endo_names,'exact'),1); //Extract net export share


//******************************simulate data from monthly model**************************************
out_withshock=NaN(M_.endo_nbr,npts,replications);//initialize matrix

//use FGRU shocks
winsorizing_dummy=1; //enable winsorizing
[shock_mat]=generate_FGRU_shocks(replications,winsorizing_dummy);

skipline();
for ii=1:replications
    //use FV et al. pruning and start at EMAS
    if ii==1
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,zeros(size(oo_.dr.ys)),zeros(M_.exo_nbr,1)); //last shock is only needed as initial condition for next run, see the description in the technical appendix of Born/Pfeifer
    else
        u_1_start=shock_mat(end,:,ii-1); //stale first order shock term from previous run used in FGRU
        [y_, y1st]= simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat(1:end-1,:,ii),options_.order,y1st(:,end),u_1_start');        
    end
    out_withshock(:,:,ii) =y_; //use all simulated time points including the leading steady state
    if mod(ii,50)==0
        fprintf('%d out of %d replications done.\n',ii,replications);
    end    
end

//Define FGRU numbers for comparison
@#if Argentina == 1 
    moments_emp(1,1)=4.77; //sigma_Y
    moments_emp(2,1)=1.31; //sigma_C/sigma_Y
    moments_emp(3,1)=3.81; //sigma_I/sigma_Y
    moments_emp(4,1)=39.33/100; //std(NX); additional 100 to counteract FV et al.'s division by 100 
    moments_emp(5,1)=-0.76; //corr(NX,Y)
    moments_emp(6,1)=moments_emp(4,1);
    moments_emp(7,1)=moments_emp(5,1);
    moments_emp(8,1)=1.78; //NX_Y at EMAS

    //Compute empirical moments
    //Data taken from FGRU empirical_moments.m
    //Sample in Argentina: 1993.Q1 - 2004.Q3, Data already seasonally adjusted

    argentina = [2016.467561	2095.869247	2137.676896	2166.951094	2175.798117	2242.036534	2239.864023	2285.238574	2246.511413	2151.566347	2136.151074	2145.778626	2232.466727	2287.22935	2337.615837	2383.442625	2441.261108	2477.999896	2550.86982	2570.334379	2591.23162	2630.790938	2607.875573	2545.079724	2500.885579	2457.383581	2440.436286	2487.918125	2499.957636	2435.335604	2424.920743	2442.370239	2445.075581	2413.789978	2301.054384	2169.392728	2043.696194	2142.344359	2131.246469	2168.339259	2216.773348	2319.914265	2360.289715	2428.737136	2483.616853	2496.180543	2560.028872
    1605.963844	1685.023459	1713.336044	1729.684839	1768.561576	1817.17525	1809.925398	1808.079891	1772.14507	1708.586507	1684.499721	1700.580324	1757.596683	1767.47444	1837.056912	1874.933674	1913.167457	1942.197719	2021.976158	2037.753747	2011.134508	2073.05066	2069.126142	2034.342373	2009.718836	2017.867655	1992.804511	2015.777775	2007.448208	1970.949082	1952.963029	1951.051781	1983.315394	1949.333817	1825.147839	1735.471017	1598.021139	1489.154227	1466.980897	1477.756631	1527.103842	1670.035294	1716.200903	1740.409184	1747.271476	1795.925268	1823.679464
    440.401957	470.1074034	498.9707866	506.4000496	500.3671204	512.1548825	510.90428	536.5915389	510.0848211	426.5299916	460.499141	462.3868639	490.7800217	525.7711258	529.1684184	531.3441593	564.4102345	596.0305228	605.1044022	615.5751767	652.28283	638.1696163	617.2442759	576.7648241	533.2269932	494.9932933	503.8938121	512.8872207	517.1609677	500.2337674	484.6050151	485.5896064	482.025197	472.3574426	413.7977432	327.912866	186.8928758	307.9665196	248.4189732	294.745173	320.2710351	374.5377063	382.5689289	455.8055897	508.3327661	508.7994571	545.9706644
    -25.53093131	-66.2630349	-68.60715319	-73.09687997	-88.05088558	-93.56859263	-75.11680056	-61.2627784	-33.8868488	8.666970808	-3.069221696	-14.06030091	-15.87221991	-17.3729521	-20.79617005	-16.81735279	-42.37149796	-67.6592032	-70.371787	-75.88774685	-82.21457327	-79.15282628	-76.33406291	-63.6779483	-50.27843905	-46.52274761	-59.62486371	-42.71847035	-27.96263341	-24.61147658	-21.62606656	2.287599415	-15.54782754	-0.025148149	54.26939716	103.9272826	247.9283239	359.1558529	412.1323747	397.6816383	354.1510176	291.8875888	259.1683081	230.8455176	210.5368151	214.8534078	186.8694046
    2020.834869	2088.867827	2143.699678	2162.988008	2180.877811	2235.76154	2245.712878	2283.408651	2248.343042	2143.783469	2141.92964	2148.906887	2232.504485	2275.872614	2345.42916	2389.460481	2435.206193	2470.569039	2556.708774	2577.441177	2581.202765	2632.06745	2610.036355	2547.429249	2492.66739	2466.3382	2437.073459	2485.946525	2496.646542	2446.571373	2415.941978	2438.928987	2449.792763	2421.666111	2293.214979	2167.311166	2032.842339	2156.276599	2127.532245	2170.183443	2201.525895	2336.46059	2357.93814	2427.060292	2466.141057	2519.578133	2556.519533] ;

    //extract cyclical properties
    [yt,yd] = hpfilter(log(argentina(1,:)'),1600) ;
    [nx_y_t,nx_y_d] = hpfilter((argentina(4,:)./(argentina(1,:)))',1600) ;
    moments_emp(9,1)=std(nx_y_d)*100;
    moments_emp(10,1)=corr(nx_y_d,yd);

    FGRU_moments(1,1)=5.30;
    FGRU_moments(2,1)=1.54;
    FGRU_moments(3,1)=3.90;
    FGRU_moments(4,1)=NaN;
    FGRU_moments(5,1)=NaN;
    FGRU_moments(6,1)=0.48;
    FGRU_moments(7,1)=0.05;
    FGRU_moments(8,1)=1.75;
    FGRU_moments(9:10,1)=NaN;
@# endif

@#if Ecuador == 1 
    moments_emp(1,1)=2.46;
    moments_emp(2,1)=2.48;
    moments_emp(3,1)=9.32;
    moments_emp(4,1)=65/100; //additional 100 to counteract FV et al.'s division by 100 
    moments_emp(5,1)=-0.60;
    moments_emp(6,1)=moments_emp(4,1);
    moments_emp(7,1)=moments_emp(5,1);
    moments_emp(8,1)=3.86;
    moments_emp(9:10,1)=NaN;

    FGRU_moments(1,1)=2.23;
    FGRU_moments(2,1)=2.13;
    FGRU_moments(3,1)=9.05;
    FGRU_moments(4,1)=NaN;
    FGRU_moments(5,1)=NaN;
    FGRU_moments(6,1)=1.77;
    FGRU_moments(7,1)=-0.04;
    FGRU_moments(8,1)=3.95;
    FGRU_moments(9:10,1)=NaN;
@# endif


@#if Venezuela == 1 
    moments_emp(1,1)=4.72;
    moments_emp(2,1)=0.87;
    moments_emp(3,1)=3.42;
    moments_emp(4,1)=18/100; //additional 100 to counteract FV et al.'s division by 100 
    moments_emp(5,1)=-0.11;
    moments_emp(6,1)=moments_emp(4,1);
    moments_emp(7,1)=moments_emp(5,1);
    moments_emp(8,1)=4.07;
    moments_emp(9:10,1)=NaN;

    FGRU_moments(1,1)=4.56;
    FGRU_moments(2,1)=0.51;
    FGRU_moments(3,1)=3.81;
    FGRU_moments(4,1)=NaN;
    FGRU_moments(5,1)=NaN;
    FGRU_moments(6,1)=1.60;
    FGRU_moments(7,1)=-0.10;
    FGRU_moments(8,1)=4.14;
    FGRU_moments(9:10,1)=NaN;
@# endif


@#if Brazil == 1 
    moments_emp(1,1)=4.64;
    moments_emp(2,1)=1.10;
    moments_emp(3,1)=1.65;
    moments_emp(4,1)=23/100; //additional 100 to counteract FV et al.'s division by 100 
    moments_emp(5,1)=-0.26;
    moments_emp(6,1)=moments_emp(4,1);
    moments_emp(7,1)=moments_emp(5,1);
    moments_emp(8,1)=0.1;
    moments_emp(9:10,1)=NaN;

    FGRU_moments(1,1)=4.52;
    FGRU_moments(2,1)=0.44;
    FGRU_moments(3,1)=1.67;
    FGRU_moments(4,1)=NaN;
    FGRU_moments(5,1)=NaN;
    FGRU_moments(6,1)=0.60;
    FGRU_moments(7,1)=0.18;
    FGRU_moments(8,1)=0.52;
    FGRU_moments(9:10,1)=NaN;
@# endif



//******************************Aggregate data and evaluate moments**************************************
//compute moments on model data
[moments_short]=get_quarterly_moments(out_withshock(:,:,1:200),ergodicmean_no_shocks,M_,oo_);
// moments_short(9,:)=[nx_y_share nx_y_share];
[moments_long, row_names, column_names, std_nx]=get_quarterly_moments(out_withshock,ergodicmean_no_shocks,M_,oo_);
// moments_long(9,:)=[nx_y_share nx_y_share];
stdevs=cumsum(std_nx)./cumsum(ones(replications,3));

//original calibration
fprintf('%25s\n',['Moments Original ',end_save_string])

fprintf('%25s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n','','Data','FGRU','Sum(200)','Sum(10000)','Mean(200)','Mean(10000)')
for ii=1:size(row_names,1)
fprintf('%25s \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f \t %10.2f\n',row_names{ii,1},moments_emp(ii,1),FGRU_moments(ii,1),moments_short(ii,1),moments_long(ii,1),moments_short(ii,2),moments_long(ii,2))
end

//******************************Generate and display IRFs**************************************

shock_mat=zeros(npts,M_.exo_nbr); //initialize shock matrix
shock_mat(1,strmatch('u_sigma_r',M_.exo_names,'exact'))=1; //set risk shock to 1 standard deviation
shock_string=deblank(M_.exo_names(strmatch('u_sigma_r',M_.exo_names,'exact'),:)); //extract name
out_withshock=simult_FGRU(ergodicmean_no_shocks,oo_.dr,shock_mat,options_.order,zeros(size(oo_.dr.ys)),shock_mat(1,:)'); //simulate series after shock
%% To simulate serires with proper pruning scheme as in Andreasen et al. (2013), use the following code
% out_withshock = simult_(ergodicmean_no_shocks,oo_.dr,shock_mat,options_.order); //simulate series
IRFS = (out_withshock - ergodicmean_no_shocks*ones(1,1+npts))*100; //express as percentage deviations from EMAS; entries that need to be treated differently are subsequently overwritten

// express interest rates as annualized basis points (12 months times 100 to be in percentage points time 100 for basis points)
IRFS(strmatch('eps_r',M_.endo_names,'exact'),:)=12*100*100*(out_withshock(strmatch('eps_r',M_.endo_names,'exact'),:)-ergodicmean_no_shocks(strmatch('eps_r',M_.endo_names,'exact'))); //yearly interest rate in basis points
IRFS(strmatch('eps_tb',M_.endo_names,'exact'),:)=12*100*100*(out_withshock(strmatch('eps_tb',M_.endo_names,'exact'),:)-ergodicmean_no_shocks(strmatch('eps_tb',M_.endo_names,'exact'))); //yearly interest rate in basis points
IRFS(strmatch('r',M_.endo_names,'exact'),:)=12*100*100*(out_withshock(strmatch('r',M_.endo_names,'exact'),:)-ergodicmean_no_shocks(strmatch('r',M_.endo_names,'exact'))); //yearly interest rate in basis points

//express debt in percentage deviations from ergodic mean in absence of shocks
IRFS(strmatch('D',M_.endo_names,'exact'),:)=IRFS(strmatch('D',M_.endo_names,'exact'),:)/(ergodicmean_no_shocks(strmatch('D',M_.endo_names,'exact'))); //deviation from debt; debt was not logged, so the difference from the EMAS needs to be normalized by EMAS


plot_vars={'C';'I';'Y';'H';'eps_r';'D'};
plot_vars_heading={'Consumption';'Investment';'Output';'Hours';'Interest Rate Spread';'Debt'};


//aggregate to quarterly IRFs
n_quarters=floor(npts/3);
for ii=1:n_quarters
    IRFS_quarterly(:,ii)=mean(IRFS(:,1+(ii-1)*3+1:1+ii*3),2); //Correct aggregation by mean; start aggregation at point 2 to get rid of initial condition
    IRFS_quarterly_FGRU(:,ii)=sum(IRFS(:,1+(ii-1)*3+1:1+ii*3),2); //FGRU aggregation by sum; start aggregation at point 2 to get rid of initial condition
end

//take care of properly aggregated variables in FGRU by dividing them by three
IRFS_quarterly_FGRU(strmatch('r',M_.endo_names,'exact'),:)=IRFS_quarterly_FGRU(strmatch('r',M_.endo_names,'exact'),:)/3;
IRFS_quarterly_FGRU(strmatch('eps_r',M_.endo_names,'exact'),:)=IRFS_quarterly_FGRU(strmatch('eps_r',M_.endo_names,'exact'),:)/3;
IRFS_quarterly_FGRU(strmatch('eps_tb',M_.endo_names,'exact'),:)=IRFS_quarterly_FGRU(strmatch('eps_tb',M_.endo_names,'exact'),:)/3;
IRFS_quarterly_FGRU(strmatch('D',M_.endo_names,'exact'),:)=IRFS_quarterly_FGRU(strmatch('D',M_.endo_names,'exact'),:)/3;


figure('name','Comparison FGRU IRFS')
for ii=1:6
    subplot(3,2,ii)
    IRF_plot=IRFS_quarterly(strmatch(plot_vars(ii,:),M_.endo_names,'exact'),:);
    IRF_plot_FGRU=IRFS_quarterly_FGRU(strmatch(plot_vars(ii,:),M_.endo_names,'exact'),:);
    if max(abs(IRF_plot)) >1e-10
        plot(1:n_quarters,IRF_plot,'b-',1:n_quarters,IRF_plot_FGRU,'r--','LineWidth',1.5)
    else
        plot(1:n_quarters,zeros(1,n_quarters),'b-',1:n_quarters,zeros(1,n_quarters),'r--','LineWidth',1.5) 
    end
    if ii==5
    legend('Correct Aggregation','FGRU Aggregation')
    end
    title(plot_vars_heading(ii,:))
    axis tight
    grid on
    hold on
    plot((1:n_quarters),zeros(n_quarters),'r')
end
eval(['print -depsc2 IRF_comparison_',country_string,end_save_string])


//******************************Generate Figure 6**************************************

//Compute IRFS, accounting for variables in percent instead of percentage points and the fact that the EMAS is 0 for CA and NX
for ii=1:n_quarters
    IRFS = (out_withshock - ergodicmean_no_shocks*ones(1,1+npts))*100;
    D_IRF_quarterly(:,ii)=mean(out_withshock(strmatch('D',M_.endo_names,'exact'),(ii-1)*3+1:ii*3));
    Y_IRF_quarterly(:,ii)=sum(exp(out_withshock(strmatch('Y',M_.endo_names,'exact'),(ii-1)*3+1:ii*3))); //get quarterly output as sum of exp() of logged monthly values
    NX_IRF_quarterly(:,ii)=sum(out_withshock(strmatch('NX',M_.endo_names,'exact'),(ii-1)*3+1:ii*3)); //absolute deviation from EMAS/SS of 0
    CA_IRF_quarterly(:,ii)=sum(out_withshock(strmatch('CA',M_.endo_names,'exact'),(ii-1)*3+1:ii*3)); //absolute deviation from EMAS/SS of 0
end

figure
//Debt to GDP ratio
subplot(3,1,1)
D_Y_quarterly_mean=ergodicmean_no_shocks(strmatch('D',M_.endo_names,'exact'),:)/(3*exp(ergodicmean_no_shocks(strmatch('Y',M_.endo_names,'exact'),:)))*100; //mean for deviation from it
D_Y_IRF=D_IRF_quarterly./Y_IRF_quarterly*100; //evolution of debt and GDP, in levels, not in deviations from SS
plot(1:n_quarters,D_Y_IRF,'b-',1:n_quarters,D_Y_quarterly_mean*ones(1,n_quarters),'r--','LineWidth',1.5)
title('D/Y','FontSize',14)
axis tight
set(gca,'FontSize',12)
//Current account to GDP ratio
subplot(3,1,2)
CA_Y_IRF_quarterly=CA_IRF_quarterly./Y_IRF_quarterly*100;
plot(1:n_quarters,CA_Y_IRF_quarterly,1:n_quarters,zeros(1,n_quarters),'r--','LineWidth',1.5)
title('CA/Y','FontSize',14)
set(gca,'FontSize',12)
axis tight
//Net Export to GDP ratio
subplot(3,1,3)
NX_Y_quarterly_mean=ergodicmean_no_shocks(strmatch('NX',M_.endo_names,'exact'),:)/exp(ergodicmean_no_shocks(strmatch('Y',M_.endo_names,'exact'),:))*100;
NX_Y_percent_IRF_quarterly=NX_IRF_quarterly./Y_IRF_quarterly*100-NX_Y_quarterly_mean;
plot(1:n_quarters,NX_Y_percent_IRF_quarterly,1:n_quarters,zeros(1,n_quarters),'r--','LineWidth',1.5)
title('NX/Y','FontSize',14)
axis tight
set(gca,'FontSize',12)
eval(['print -depsc2 Current_Account',country_string,end_save_string])


@#if estimation == 1 
    x_start=[sigma_x,phipar,D_bar,exp(Phi)]; //use calibration as starting point
    target=[moments_emp(1:3);moments_emp(8)]'; //define target moments
    //optimizer options
    H0 = 1e-2*eye(length(x_start)); //Initial Hessian 
    crit = 1e-7; //Tolerance
    nit = 1000; //Number of iterations
    //make sure Dynare does not print out stuff during runs
    options_.nocorr=1;
    options_.noprint=1;
    options_.verbosity=0;

    //options_.qz_criterium = 1+1e-6; //required because it is empty by default, leading to a crash in k_order_pert
    [shock_mat]=generate_FGRU_shocks(estimation_replications,winsorizing_dummy);
    [fhat,xhat] = csminwel(@smm_diff_function,x_start,H0,[],crit,nit,target,estimation_replications,shock_mat(:,:,1:estimation_replications));

@# endif

clear  out_withshock shock_mat out_noshock winsorizing temp log_deviations_SS_noshock; //delete big matrices before saving
save([country_string,'_',end_save_string]); //save results
