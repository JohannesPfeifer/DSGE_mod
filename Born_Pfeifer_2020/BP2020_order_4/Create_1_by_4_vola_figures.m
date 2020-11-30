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

h_tfp_vola=figure('Name','Wedge Responses, Technology Shock Volatility');
set(h_tfp_vola,'PaperType','a4','PaperPositionMode','manual','PaperUnits','centimeters','PaperPosition',[1,1,28,21],'renderer', 'painters');
subplot(4,4,1)
plot(1:options_.irf,100*oo_.irfs.sigma_z_eps_sigma_z,'LineWidth',1.5)
title('TFP Uncertainty','FontSize',10)
xlabel('quarters','FontSize',10)
ylabel('percent','FontSize',10)
set(gca,'FontSize',10)
axis tight
hline(0);

subplot(4,4,2)
if max(abs(oo_.irfs.firm_wedge_eps_sigma_z))<1e-12
    oo_.irfs.firm_wedge_eps_sigma_z=zeros(size(oo_.irfs.firm_wedge_eps_sigma_z));
end  
plot(1:options_.irf,100*oo_.irfs.firm_wedge_eps_sigma_z,'LineWidth',1.5)
title('Price Markup','FontSize',10)
xlabel('quarters','FontSize',10)
set(gca,'FontSize',10)
axis tight
hline(0);

subplot(4,4,3)
if max(abs(oo_.irfs.household_wedge_eps_sigma_z))<1e-12
    oo_.irfs.household_wedge_eps_sigma_z=zeros(size(oo_.irfs.household_wedge_eps_sigma_z));
end  
plot(1:options_.irf,100*oo_.irfs.household_wedge_eps_sigma_z,'LineWidth',1.5)
title('Wage Markup','FontSize',10)
xlabel('quarters','FontSize',10)
set(gca,'FontSize',10)
axis tight
hline(0);

subplot(4,4,4)
plot(1:options_.irf,100*oo_.irfs.log_Y_eps_sigma_z,'LineWidth',1.5)
title('Output','FontSize',10)
xlabel('quarters','FontSize',10)
set(gca,'FontSize',10)
axis tight
hline(0);

h_G_vola=figure('Name','Wedge Responses, Government Spending Shock Volatility');
set(h_G_vola,'PaperType','a4','PaperPositionMode','manual','PaperUnits','centimeters','PaperPosition',[1,1,28,21],'renderer', 'painters');
subplot(4,4,1)
plot(1:options_.irf,100*oo_.irfs.sigma_G_eps_sigma_G,'LineWidth',1.5)
title('G Uncertainty','FontSize',10)
axis tight
xlabel('quarters','FontSize',10)
ylabel('percent','FontSize',10)
set(gca,'FontSize',10)
hline(0);

subplot(4,4,2)
if max(abs(oo_.irfs.firm_wedge_eps_sigma_G))<1e-12
    oo_.irfs.firm_wedge_eps_sigma_G=zeros(size(oo_.irfs.firm_wedge_eps_sigma_G));
end  
plot(1:options_.irf,100*oo_.irfs.firm_wedge_eps_sigma_G,'LineWidth',1.5)
title('Price Markup','FontSize',10)
xlabel('quarters','FontSize',10)
set(gca,'FontSize',10)
axis tight
hline(0);

subplot(4,4,3)
if max(abs(oo_.irfs.household_wedge_eps_sigma_G))<1e-12
    oo_.irfs.household_wedge_eps_sigma_G=zeros(size(oo_.irfs.household_wedge_eps_sigma_G));
end  
plot(1:options_.irf,100*oo_.irfs.household_wedge_eps_sigma_G,'LineWidth',1.5)
title('Wage Markup','FontSize',10)
xlabel('quarters','FontSize',10)
set(gca,'FontSize',10)
axis tight
hline(0);

subplot(4,4,4)
plot(1:options_.irf,100*oo_.irfs.log_Y_eps_sigma_G,'LineWidth',1.5)
title('Output','FontSize',10)
xlabel('quarters','FontSize',10)
set(gca,'FontSize',10)
axis tight
hline(0);
