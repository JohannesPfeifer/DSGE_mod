load gpr.dat;
timeline=(1948.25:0.25:2003)';
%% Post 1980s
gobs = demean(gpr(timeline>=1980,1));
piobs = demean(gpr(timeline>=1980,2));
robs = demean(gpr(timeline>=1980,3));