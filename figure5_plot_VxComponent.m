
clear;clc;close all
load('figure5_conventional_13ms.mat')
plotimage(-Seismic_Vx(:,45:end-45));
title('')
xlabel('x/dx')
ylabel('travel time/dt')



load('figure5_HEI_13ms.mat')
plotimage(-Seismic_Vx(:,45:end-45));
title('')
xlabel('x/dx')
ylabel('travel time/dt')

load('figure5_conventional_26ms.mat')
plotimage(-Seismic_Vx(:,45:end-45));
title('')
xlabel('x/dx')
ylabel('travel time/dt')

load('figure5_HEI_26ms.mat')
plotimage(-Seismic_Vx(:,45:end-45));
title('')
xlabel('x/dx')
ylabel('travel time/dt')