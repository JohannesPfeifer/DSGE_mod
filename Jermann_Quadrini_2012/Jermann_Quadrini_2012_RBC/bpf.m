function yf=bpf(y,up,dn,K)
% function yf=bpf(y,up,dn,K)
% Program to compute band-pass filtered series
% Inputs are
%  y:   data (rows = observations, columns=series)
%  up:  period corresponding to highest frequency (e.g., 6)
%  dn:  period corresponding to lowest frequency (e.g., 32)
%  K:   number of terms in approximating moving average
%  [calls filtk.m (filter with symmetric weights) as subroutine]
%
% based on the baxter.m at http://users.cla.umn.edu/~erm/data/sr338/replicate/baxter.m
%
% Copyright: Ellen MacGrattan 2004

x=[up dn];

if (up>dn)
    disp('Periods reversed: switching indices up & dn')
    disp(' ')
    dn=x(1); up=x(2);
end

if (up<2)
    up=2;
    disp('Higher periodicity > max: Setting up=2')
    disp(' ')
end

% convert to column vector
[r c]=size(y);
if (r<c)
    y=y';
    disp('There are more columns than rows: Transposing data matrix')
    disp(' ')
end

% Implied Frequencies
omubar=2*pi/up;
omlbar=2*pi/dn;

% An approximate low pass filter, with a cutoff frequency of "ombar",
% has a frequency response function
%
%  alpha(om) = a0 + 2*a1 cos(om) + ... 2*aK cos(K om)
%
% and the ak's are given by:
%
% a0 = ombar/(pi)             ak = sin(k ombar)/(k pi)
%
% where ombar is the cutoff frequency.

% A band-pass filter is the difference between two
% low-pass filters,
%   bp(L) = bu(L) - bl(L)
% with bu(L) being the filter with the high cutoff point and bl(L) being
% that with the low cutoff point.  Thus, the weights are differences
% of weights for two low-pass filters.

% Construct filter weights for bandpass filter (a(0)....a(K)).

akvec=NaN(1,K+1);

akvec(1)=(omubar-omlbar)/(pi);  % weight at k=0

for k=1:K;
    akvec(k+1)=(sin(k*omubar)-sin(k*omlbar))/(k*pi); % weights at k=1,2,...K
end

% Impose constraint on frequency response at om = 0
% (If high pass filter, this amounts to requiring that weights sum to zero).
% (If low pass filter, this amounts to requiring that weights sum to one).

if (dn>1000)
    disp('dn > 1000: assuming low pass filter')
    phi=1;
else
    phi=0;
end

% sum of weights without constraint
theta=akvec(1)+2*sum(akvec(2:K+1));
% amount to add to each nonzero lag/lead to get sum = phi
theta=(phi-theta)/(2*K+1);
% adjustment of weights
akvec=akvec+theta;

% filter the time series
yf=filtk(y,akvec);

if (r<c)
    yf=yf';
end
end

function  yf=filtk(y,a)

% Filter data with a filter with symmetric filter with weights
% data is organized (rows=obs, columns=series)
% a=[a0 a1 ... aK];

K=max(size(a))-1;  % max lag;

T=max(size(y));    % number of observations;

% Set vector of weights

avec=NaN(1,2*K+1);
avec(K+1)=a(1);
for i=1:K;
    avec(K+1-i)=a(i+1);
    avec(K+1+i)=a(i+1);
end
yf=zeros(size(y));
for t=K+1:1:T-K
    yf(t,:)=avec*y(t-K:t+K,:);
end
yf(1:K,:)=NaN;
yf(T-K+1:T,:)=NaN;

end