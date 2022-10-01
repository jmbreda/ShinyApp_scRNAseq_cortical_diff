clear all; clc;

addpath('~/Matlab_Scripts/');

cluster_colors = distinguishable_colors(100);
[fid,msg] = fopen('../Zahra_dev/data/cluster_colors.txt','w')
for i = 1:size(cluster_colors,1);
	fprintf(fid,'%s\n',rgb2hex(cluster_colors(i,:)));
end
fclose(fid);

time_colors = [[1 0 0];[1 .6 0];[.13 .55 .13];[0 .8 1];[0 0 1]];
[fid,msg] = fopen('../Zahra_dev/data/time_colors.txt','w')
for i = 1:size(time_colors,1);
	fprintf(fid,'%s\n',rgb2hex(time_colors(i,:)));
end
fclose(fid);
