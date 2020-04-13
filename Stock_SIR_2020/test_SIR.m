gamma=0.55;
pi_a=0.3;
s_0=0.02;
N_total=329000000;
T=52;
S=NaN(T+1,1);
R=NaN(T+1,1);
I=NaN(T+1,1);
P_I_and_Symptoms=NaN(T+1,1);
P_symptomatic=NaN(T+1,1);
P_I_given_symptoms=NaN(T+1,1);
Total_infected=NaN(T+1,1);
R_0=NaN(T+1,1);
N=NaN(T+1,1);

beta_temp=NaN(T,1);
beta_temp(1:12)=2.1;
beta_temp(13:24)=0.66;
beta_temp(25:28)=0.8800;
beta_temp(29:32)=0.9900;
beta_temp(33:36)=1.1;
beta_temp(37:40)=1.21;
beta_temp(41:44)=1.32;
beta_temp(45:T)=2.1;
beta=NaN(T+1,1);
beta(2:end)=beta_temp;

I(1)=50/N_total;
R(1)=0/N_total;
S(1)=1-I(1)-R(1);

for ii=2:T+1
    S(ii)=S(ii-1)-beta(ii)*I(ii-1)*S(ii-1);
    R(ii)=R(ii-1)+gamma*I(ii-1);
    I(ii)=I(ii-1)+beta(ii)*I(ii-1)*S(ii-1)-gamma*I(ii-1);
    P_I_and_Symptoms(ii)=(1-pi_a)*I(ii);
    P_symptomatic(ii)=P_I_and_Symptoms(ii)+s_0*(S(ii)+R(ii));
    P_I_given_symptoms(ii)=(1-pi_a)*I(ii)/P_symptomatic(ii);
    Total_infected(ii)=I(ii)+R(ii);
    R_0(ii)=beta(ii)/gamma;
    N(ii)=S(ii)+I(ii)+R(ii);
end