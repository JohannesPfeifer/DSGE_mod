function plot_NBER_recessions(plothandle,legend_strings)
% function plot_NBER_recessions(plothandle,legend_strings)
% Plots NBER recessions as shaded areas into background of a graph
% Inputs:
%   plothandle:         [array of plot-handles]     existing plot handles
%                                                   that should be in foreground
%   legend_strings      [character array]           optional legend to
%                                                   print
% Outputs: None

%
% Copyright (C) 2013-2016 Benjamin Born and Johannes Pfeifer
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

recessiondates=[1953.25,1954.25;
    1957.50,1958.25;
    1960.25,1961;
    1969.75,1970.75;
    1973.75,1975;
    1980,1980.50;
    1981.50,1982.75;
    1990.50,1991;
    2001,2001.75;
    2007.75,2009.25];

axis tight
aaa=ylim;
bottom=aaa(1,1);
top=aaa(1,2);

%delete dates that happen before first data point
x_dates=xlim;
recessiondates(recessiondates(:,2)<x_dates(1),:)=[];

for ii=1:length(recessiondates)
    hold on
    ha = area([recessiondates(ii,1) recessiondates(ii,2)], [bottom top-bottom; bottom top-bottom],'FaceColor',[0.9 0.9 0.9],'EdgeColor','white','ShowBaseline','off');
    set(ha(1), 'FaceColor', 'none') % this makes the bottom area invisible
    set(ha, 'LineStyle', '-')
    set(ha, 'LineStyle', '-')
end
ylim(aaa)

if nargin==3
    legend(plothandle,legend_strings);
end
hline(0,'k-')
uistack(plothandle,'top')
set(gca,'Layer','top')


function hline(y,linetype)

g=ishold(gca);
hold on

x=get(gca,'xlim');
h=plot(x,[y y],linetype);
if g==0
    hold off
end
set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends