function dyn_latex_table_modified(options_,title,LaTeXtitle,headers,labels,values,label_width,val_width,val_precis,optional_header)
%function dyn_latex_table_modified(options_,title,LaTeXtitle,headers,labels,values,label_width,val_width,val_precis,optional_header)

% Copyright (C) 2015-2017 Dynare Team
% Copyright (C) 2018 Johannes Pfeifer (removed code that prevents calling function from base workspace)
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if options_.noprint
    return
end

%% get width of label column
if ~isempty(label_width)
    label_width = max(size(deblank(char(headers(1,:),labels)),2)+2, ...
                      label_width);
else %use default length
    label_width = max(size(deblank(char(headers(1,:),labels)),2))+2;
end
label_format_leftbound  = sprintf('$%%-%ds$',label_width);

%% get width of label column
if all(~isfinite(values))
    values_length = 4;
else
    values_length = max(ceil(max(max(log10(abs(values(isfinite(values))))))),1)+val_precis+1;
end
if any(values) < 0 %add one character for minus sign
    values_length = values_length+1;
end

%% get width of header strings
headers_length = max(size(deblank(headers(2:end,:)),2));
if ~isempty(val_width)
    val_width = max(max(headers_length,values_length)+4,val_width);
else
    val_width = max(headers_length,values_length)+4;
end
value_format  = sprintf('%%%d.%df',val_width,val_precis);
header_string_format  = sprintf('$%%%ds$',val_width);

%Create and print header string
if length(headers) > 0
    header_string = sprintf(label_format_leftbound ,strrep(deblank(headers(1,:)),'\','\\'));
    header_code_string='l';
    for i=2:size(headers,1)
        header_string  = [header_string '\t & \t ' sprintf(header_string_format,strrep(deblank(headers(i,:)),'\','\\'))];
        header_code_string= [header_code_string 'c'];
    end
end
header_string=[header_string '\\\\\n'];

filename = [LaTeXtitle '.tex'];
fidTeX = fopen(filename,'w');

stack = dbstack;
% fprintf(fidTeX,['%% ' datestr(now,0) ', created by ' stack(2).file]);
fprintf(fidTeX,' \n');
fprintf(fidTeX,' \n');
fprintf(fidTeX,'\\begin{center}\n');
fprintf(fidTeX,['\\begin{longtable}{%s} \n'],header_code_string);
fprintf(fidTeX,['\\caption{',title,'}\\\\\n ']);

fprintf(fidTeX,['\\label{Table:',LaTeXtitle,'}\\\\\n']);
fprintf(fidTeX,'\\toprule \n');
if nargin==11
    for ii=1:size(optional_header,1)
        fprintf(fidTeX,'%s\n',optional_header{ii});
    end
end
fprintf(fidTeX,header_string);
fprintf(fidTeX,'\\midrule \\endfirsthead \n');
fprintf(fidTeX,'\\caption{(continued)}\\\\\n ');
fprintf(fidTeX,'\\toprule \\\\ \n');
if nargin==11
    for ii=1:size(optional_header,1)
        fprintf(fidTeX,'%s\n',optional_header{ii});
    end
end
fprintf(fidTeX,header_string);
fprintf(fidTeX,'\\midrule \\endhead \n');
fprintf(fidTeX,['\\midrule \\multicolumn{',num2str(size(headers,1)),'}{r}{(Continued on next page)} \\\\ \\bottomrule \\endfoot \n']);
fprintf(fidTeX,'\\bottomrule \\endlastfoot \n');
for i=1:size(values,1)
    fprintf(fidTeX,label_format_leftbound,deblank(labels(i,:)));
    fprintf(fidTeX,['\t & \t' value_format],values(i,:));
    fprintf(fidTeX,' \\\\ \n');
end

fprintf(fidTeX,'\\end{longtable}\n ');
fprintf(fidTeX,'\\end{center}\n');
fprintf(fidTeX,'%% End of TeX file.\n');
fclose(fidTeX);