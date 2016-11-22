function C = nancorr( X,Y )
% NANCORR calculates the sample correlation coefficient
%    for the series with NaNs expected.
%    X is the one series, Y is another.

% Source: https://www.mathworks.com/matlabcentral/fileexchange/35501-nancorr/content/nancorr.m
% Copyright (c) 2012, Isabella Osetinsky
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

X=X(:);
Y=Y(:);
L1=length(X);
L2=length(Y);

if L1 ~= L2
    error('The samples must be of the same length')
end

for i=1:L1,
    if isnan(X(i)),
        Y(i)=nan;
    end
    if isnan(Y(i)),
        X(i)=nan;
    end
end
        
Xm=nanmean(X);
Ym=nanmean(Y);
C=nansum((X-Xm).*(Y-Ym))/sqrt((nansum((X-Xm).^2))*(nansum((Y-Ym).^2)));
end