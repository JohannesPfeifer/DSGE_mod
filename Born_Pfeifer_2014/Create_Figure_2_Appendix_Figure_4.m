% This file generates Figure 2. Convergence Behavior of the Relative Net Export Volatility Statistic
% and Appendix Figure 4. Convergence Behavior of Different Net Export Volatility Statistics in the Re-
% calibrated Model of:
% Benjamin Born and Johannes Pfeifer (2014): "Risk Matters: A comment", American Economic Review
% 
% Requirements: requires that the main mod-file has been run before in
% replications and calibration mode and that the results have been
% correctly saved. See the ReadMe-file.
% 
% Copyright (C) 2013-14 Benjamin Born and Johannes Pfeifer
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
% For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.

clear all
%specify for which calibration to plot the figure 
% end_save_string='replication'; %'replication' for Figure 1 
end_save_string='recalibration'; %'recalibration' for Appendix Figure 4 

%% load the data for Argentina
x=load(['Argentina_',end_save_string],'std_nx','moments_long','moments_short','moments_emp');

%% compute the average standard deviation up to the current iteration
NX_vol=cumsum(x.std_nx(:,1))./cumsum(ones(size(x.std_nx(:,1))));
Y_vol=cumsum(x.std_nx(:,3))./cumsum(ones(size(x.std_nx(:,3))));
stdevs_NX_over_Y=NX_vol./Y_vol; %average relative volatility

figure
subplot(2,1,1)
x_ax=1:length(stdevs_NX_over_Y);
plot(x_ax,stdevs_NX_over_Y,'b',x_ax,ones(size(x_ax))*x.moments_emp(6),'k-.','LineWidth',1.5)
text(7000,1.1*x.moments_emp(6),['Data: ', num2str(x.moments_emp(6),'%3.2f')],'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','bottom')
title('\sigma_{NX}/\sigma_Y','FontSize',12)
ylabel('Average over Rep.','FontSize',12)
axis tight
y_limits=ylim;
ylim([0 1.05*y_limits(2)])
if strcmp(end_save_string,'replication')
annotation('textarrow',[0.190191387559809 0.145933014354067],...
    [0.665529010238908 0.691126279863482],'TextEdgeColor','none',... 
    'String',num2str(stdevs_NX_over_Y(200),'%3.2f'),'FontSize',12);
elseif strcmp(end_save_string,'recalibration')
    annotation('textarrow',[0.190191387559809 0.145933014354067],...
    [0.701365187713311 0.769624573378841],'TextEdgeColor','none',... 
        'String',num2str(stdevs_NX_over_Y(200),'%3.2f'),'FontSize',12);
    
end

%% plot the convergence behavior the standard deviation of the net-export to out ratio, stored in the second column of std_nx
if strcmp(end_save_string,'recalibration')
    subplot(2,1,2)
    stdevs_NX_Y=cumsum(x.std_nx(:,2))./cumsum(ones(size(x.std_nx(:,2))))*100;
    plot(x_ax,stdevs_NX_Y,'b',x_ax,ones(size(x_ax))*x.moments_emp(10),'k-.','LineWidth',1.5)
    text(7000,x.moments_emp(10),['Data: ', num2str(x.moments_emp(10),'%3.2f')],'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','bottom')
    title('\sigma_{NX/Y}','FontSize',12)
    axis tight
    y_limits=ylim;
    ylim([0*y_limits(1) 1.05*y_limits(2)])
    if strcmp(end_save_string,'replication')
    annotation('textarrow',[0.190191387559809 0.144736842105263],...
        [0.26962457337884 0.317406143344711],'TextEdgeColor','none',... 
        'String',num2str(stdevs_NX_Y(200),'%5.2f'),'FontSize',12);
    elseif strcmp(end_save_string,'recalibration')
    annotation('textarrow',[0.190191387559809 0.144736842105263],...
        [0.284982935153584 0.332764505119455],'TextEdgeColor','none',... 
        'String',num2str(stdevs_NX_Y(200),'%5.2f'),'FontSize',12);
    end
    xlabel('Repetitions','FontSize',12)
    ylabel('Average over Rep.','FontSize',12)
end
if strcmp(end_save_string,'replication')
    print('-depsc2',['Figure_2_NX_Y_volatility_',end_save_string]);
elseif strcmp(end_save_string,'recalibration')
    print('-depsc2',['Figure_4_Appendix_NX_Y_volatility_',end_save_string]);
end
    