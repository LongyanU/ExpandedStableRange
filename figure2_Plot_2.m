
clear;clc
close all;

 load('figure2_2ms_conventional.mat')
plotimage(-ptemp3)
xlabel('x/dx')
ylabel('z/dz')
title('')

load('figure2_4ms_conventional.mat')
plotimage(-ptemp3)
xlabel('x/dx')
ylabel('z/dz')
title('')

load('figure2_2ms_HEI.mat')
plotimage(-ptemp3)
xlabel('x/dx')
ylabel('z/dz')
title('')


load('figure2_4ms_HEI.mat')

plotimage(-ptemp3)
xlabel('x/dx')
ylabel('z/dz')
title('')

