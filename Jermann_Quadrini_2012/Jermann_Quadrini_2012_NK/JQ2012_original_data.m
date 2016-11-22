% This function reads in the original Jermann/Quadrini (2012)-data and transforms it into
% demeaned growth rates used for model estimation

% Copyright (C) 2016 Johannes Pfeifer
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

startdate=1984;
enddate=2010.25;
[numbers,text]=xlsread('JQ2012.xlsx','4-Data For Estimation','A3:I257');
timeline=numbers(:,1);
if timeline(1,1)~=1947 && timeline(end,1)~=2010.25
    error('Wrong dimensions')
end

y_obs=diff(log(numbers(:,strcmp('Gross domestic product',text(1,:)))));
c_obs=diff(log(numbers(:,strcmp('Personal consumption expenditures',text(1,:)))));
invest_obs=diff(log(numbers(:,strcmp('Gross private domestic investment',text(1,:)))));
pi_obs=diff(log(numbers(:,strcmp('Implicit price deflator for GDP',text(1,:)))));
r_obs=log(numbers(2:end,strcmp('Federal fund rate (gross, quarterly rate)',text(1,:))));
n_obs =diff(log(numbers(:,strcmp('Hours',text(1,:)))));
W_obs=diff(log(numbers(:,strcmp('Wages',text(1,:)))));
debt_repurchase_obs=numbers(2:end,strcmp('Debt Repurchase',text(1,:)));
timeline=timeline(2:end);

y_obs=demean(y_obs(timeline>=startdate & timeline<=enddate));
c_obs=demean(c_obs(timeline>=startdate & timeline<=enddate));
invest_obs=demean(invest_obs(timeline>=startdate & timeline<=enddate));
pi_obs=demean(pi_obs(timeline>=startdate & timeline<=enddate));
r_obs=demean(r_obs(timeline>=startdate & timeline<=enddate));
n_obs =demean(n_obs(timeline>=startdate & timeline<=enddate));
W_obs=demean(W_obs(timeline>=startdate & timeline<=enddate));
debt_repurchase_obs=demean(debt_repurchase_obs(timeline>=startdate & timeline<=enddate));
timeline=timeline(timeline>=startdate & timeline<=enddate);

% figure
% plot(timeline,y_obs)
% title('yobs')
% figure
% plot(timeline,c_obs)
% title('cobs')
% figure
% plot(timeline,invest_obs)
% title('zetainvest_obs')
% figure
% plot(timeline,pi_obs)
% title('piobs')
% figure
% plot(timeline,r_obs)
% title('robs')
% figure
% plot(timeline,n_obs)
% title('nobs')
% figure
% plot(timeline,W_obs)
% title('W_obs')
% figure
% plot(timeline,debt_repurchase_obs)
% title('debt_repurchase_obs')
