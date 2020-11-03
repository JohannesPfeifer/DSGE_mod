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

%% Do Isoelastic utility with habits
dynare BP2020_CES -Dsticky_prices=1 -Dsticky_wages=1 -DIso_elastic_utility=1
saveas(h_tfp_vola,'Figures/BP2020_tfp_vola_Iso_habits')
saveas(h_G_vola,'Figures/BP2020_G_vola_Iso_habits')

%% Do Isoelastic utility without habits
set_param_value('h',0);
        
info = stoch_simul(var_list_);
if info>0
    error('Model could not be solved')
else
    Create_1_by_4_vola_figures
end

saveas(h_tfp_vola,'Figures/BP2020_tfp_vola_Iso')
saveas(h_G_vola,'Figures/BP2020_G_vola_Iso')


%% Do Isoelastic utility without SW (2007) calibration
set_param_value('h',0.71);
set_param_value('siggma',1.4);
set_param_value('Frisch_target',0.5);
        
info = stoch_simul(var_list_);
if info>0
    error('Model could not be solved')
else
    Create_1_by_4_vola_figures
end

saveas(h_tfp_vola,'Figures/BP2020_tfp_vola_Iso_SW')
saveas(h_G_vola,'Figures/BP2020_G_vola_Iso_SW')


figure_names={'BP2020_tfp_vola_sticky_both','BP2020_tfp_vola_Iso','BP2020_tfp_vola_Iso_habits','BP2020_tfp_vola_Iso_SW'};
line_styles={'-','--',':','-.'};
f = figure('Name','TFP Vola Preferences');
set(f,'PaperType','a4','PaperPositionMode','manual','PaperUnits','centimeters','PaperPosition',[1,1,28,21],'renderer', 'painters');


for fig_iter=1:length(figure_names)
    h=openfig(['Figures/',figure_names{fig_iter}]);
    for panel_iter=1:3
        figure(f)
        target_handle = subplot(3,3,panel_iter);
        axes(target_handle);
        hold(target_handle,'on');
        pp=plot(h.Children(3-(panel_iter-1)).Children.XData,h.Children(3-(panel_iter-1)).Children.YData,'LineWidth',1.5);
        pp.LineStyle=line_styles{fig_iter};
        switch panel_iter
            case 1
                target_handle.Title.String='Price markup';
            case 2
                target_handle.Title.String='Wage markup';
            case 3
                target_handle.Title.String='Output';                
       end
        if fig_iter==length(figure_names) 
            if panel_iter==3
                legend({'Baseline','Iso','Iso Habits','SW (2007)'},'AutoUpdate','off','box','off','Location','SouthEast')
                ylim([-0.004 0.0002])
                yticks('auto')
            end
            xlim([1 20])
            hline(0);
            xlabel('quarters')
            ylabel('percent')
            box(target_handle,'on');
        end
    end
end
f.Children(2).YRuler.Exponent=-2;

saveas(f,'Figures/BP2020_tfp_vola_preferences')
print(f,'Figures/BP2020_tfp_vola_preferences','-depsc2')


figure_names={'BP2020_G_vola_sticky_both','BP2020_G_vola_Iso','BP2020_G_vola_Iso_habits','BP2020_G_vola_Iso_SW'};
line_styles={'-','--',':','-.'};
f = figure('Name','G Vola Preferences');
set(f,'PaperType','a4','PaperPositionMode','manual','PaperUnits','centimeters','PaperPosition',[1,1,28,21],'renderer', 'painters');


for fig_iter=1:length(figure_names)
    h=openfig(['Figures/',figure_names{fig_iter}]);
    for panel_iter=1:3
        figure(f)
        target_handle = subplot(3,3,panel_iter);
        axes(target_handle);
        hold(target_handle,'on');
        pp=plot(h.Children(3-(panel_iter-1)).Children.XData,h.Children(3-(panel_iter-1)).Children.YData,'LineWidth',1.5);
        pp.LineStyle=line_styles{fig_iter};
        switch panel_iter
            case 1
                target_handle.Title.String='Price markup';
            case 2
                target_handle.Title.String='Wage markup';
            case 3
                target_handle.Title.String='Output';                
       end
        if fig_iter==length(figure_names) 
            if panel_iter==3
                legend({'Baseline','Iso','Iso Habits','SW (2007)'},'AutoUpdate','off','box','off','Location','SouthEast')
                ylim([-0.00011 0.00002])
                yticks('auto')
            end
            xlim([1 20])
            hline(0);
            xlabel('quarters')
            ylabel('percent')
            box(target_handle,'on');
        end
    end
end
f.Children(2).YRuler.Exponent=-3;

saveas(f,'Figures/BP2020_G_vola_preferences')
print(f,'Figures/BP2020_G_vola_preferences','-depsc2')
