% This file creates "Figure 1: Impulse response functions to 1 percentage point 
% (annualized) monetary policy shock under Calvo" of Born/Pfeifer (2018): 
% "The New Keynesian Wage Phillips Curve: Calvo vs. Rotemberg". 
% It requires Dynare 4.5 to be in the path.

% Copyright (C) 2018 Johannes Pfeifer and Benjamin Born
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% For a copy of the GNU General Public License,
% see <http://www.gnu.org/licenses/>.

dynare Born_Pfeifer_2018_MP -DSGU_framework=0 -DCalvo=1 -Dfixed_WPC_slope=0
oo_EHL=oo_;
dynare Born_Pfeifer_2018_MP -DSGU_framework=1 -DCalvo=1 -Dfixed_WPC_slope=0
oo_SGU=oo_;

hh=figure('Name','Dynamic Responses to monetary policy shock');
subplot(2,2,1)
plot(1:options_.irf,oo_EHL.irfs.y_eps_nu,'-',1:options_.irf,oo_SGU.irfs.y_eps_nu,'--','LineWidth',2)
xlim([1 options_.irf]);
ylabel('percent')
title('Output')
ll=legend('EHL','SGU');
set(ll,'Location','SouthEast');
subplot(2,2,2)
plot(1:options_.irf,oo_EHL.irfs.pi_p_ann_eps_nu,'-',1:options_.irf,oo_SGU.irfs.pi_p_ann_eps_nu,'--','LineWidth',2)
xlim([1 options_.irf]);
title('Price inflation')
subplot(2,2,3)
plot(1:options_.irf,oo_EHL.irfs.pi_w_ann_eps_nu,'-',1:options_.irf,oo_SGU.irfs.pi_w_ann_eps_nu,'--','LineWidth',2)
xlim([1 options_.irf]);
xlabel('quarter')
ylabel('percent')
title('Wage inflation')
subplot(2,2,4)
plot(1:options_.irf,oo_EHL.irfs.w_real_eps_nu,'-',1:options_.irf,oo_SGU.irfs.w_real_eps_nu,'--','LineWidth',2)
xlim([1 options_.irf]);
xlabel('quarter')
title('Real wage')
set(findall(hh,'-property','FontWeight'),'FontWeight','normal')
print -depsc2 Born_Pfeifer_2018_MP_IRF