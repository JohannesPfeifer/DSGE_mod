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

figure_names={'BP2020_tfp_vola_sticky_both','BP2020_tfp_vola_Firm_RA','BP2020_tfp_vola_FV','BP2020_tfp_vola_FV_LL'};
f = figure('Name','TFP Vola Size effect');

figure(f)
  
target_handle = subplot(2,2,1);    
for fig_iter=1:3
    h=openfig(['Figures/',figure_names{fig_iter}]);
    axes(target_handle);
    hold(target_handle,'on')
    pp=plot(h.Children(1).Children.XData,h.Children(1).Children.YData,'LineWidth',1.5);
    if fig_iter==2
        pp.LineStyle='--';
    elseif fig_iter==3
        pp.LineStyle=':';
    end
end
f.Children(1).YRuler.Exponent=-2;
target_handle.Title.String='Output';
legend({'Baseline','Firm RA','RA + FV et al. calib.'},'AutoUpdate','off','box','off','Location','SouthEast')
xlim([1 20])
ylim([-0.01 0.002])
yticks('auto')
hline(0);
xlabel('quarters')
ylabel('percent')
box(target_handle,'on');

target_handle = subplot(2,2,2);    
for fig_iter=4
    h=openfig(['Figures/',figure_names{fig_iter}]);
    axes(target_handle);
    hold(target_handle,'on')
    pp=plot(h.Children(1).Children.XData,h.Children(1).Children.YData,'LineWidth',1.5);
    pp.Color=[0.4660    0.6740    0.1880];
    pp.LineStyle='-.';
end
f.Children(1).YRuler.Exponent=0;
target_handle.Title.String='Output';
legend({'RA + FV et al. + Leduc Liu'},'AutoUpdate','off','box','off','Location','SouthEast')
xlim([1 20])
ylim([-0.20 0.02])
yticks('auto')
hline(0);
xlabel('quarters')
ylabel('percent')
box(target_handle,'on');
  

set(f,'PaperType','a4','PaperPositionMode','manual','PaperUnits','centimeters','PaperPosition',[1,1,28,21],'renderer', 'painters');
set(findall(f,'-property','ShowBaseLine'),'ShowBaseLine','Off');
set(findall(f,'-property','FontSize'),'FontSize',12)
set(findall(f,'-property','FontWeight'),'FontWeight','normal')

print(f,'Figures/BP2020_tfp_vola_effect_size_order_4','-depsc2')



f = figure('Name','TFP Vola Size effect order 4');

figure(f)
  
target_handle = subplot(1,1,1);    
h=openfig(['Figures/BP2020_tfp_vola_FV_LL']);
axes(target_handle);
hold(target_handle,'on')
% pp=plot(h.Children(1).Children.XData,h.Children(1).Children.YData,'LineWidth',1.5);
pp=plot(h.Children(1).Children.XData,2*h.Children(1).Children.YData,'LineWidth',1.5);
h=openfig(['Figures/BP2020_tfp_vola_FV_LL_large']);
axes(target_handle);
hold(target_handle,'on')
pp=plot(h.Children(1).Children.XData,h.Children(1).Children.YData,'LineWidth',1.5);
pp.LineStyle='--';
h=openfig(['Figures/BP2020_tfp_vola_FV_LL_large_sequential']);
axes(target_handle);
hold(target_handle,'on')
pp=plot(h.Children(1).Children.XData,h.Children(1).Children.YData,'LineWidth',1.5);
pp.LineStyle=':';
f.Children(1).YRuler.Exponent=0;
target_handle.Title.String='Output';
legend({'RA + FV et al. + Leduc Liu\newline 2\sigma times 2','RA + FV et al. + Leduc Liu\newline 4\sigma','RA + FV et al. + Leduc Liu\newline 2\sigma at t=1,2'},'AutoUpdate','off','box','off','Location','SouthEast')
xlim([1 20])
ylim([-0.45 0.1])
yticks('auto')
hline(0);
xlabel('quarters')
ylabel('percent')
box(target_handle,'on');

print('Figures/BP2020_tfp_vola_effect_size_order_4_comp','-depsc2')