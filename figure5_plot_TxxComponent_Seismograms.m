
clear;clc;close all

load('figure5_conventional_13ms.mat')
figure;plot(-Seismic_Txx(:,250),'r','linewidth', 1);
temp=Seismic_Txx;
grid on
axis([0 4000 -12 7])

load('figure5_HEI_13ms.mat')
hold on;plot(-Seismic_Txx(:,250),'k','linewidth', 1);
hold on;plot((-Seismic_Txx(:,250)+temp(:,250)),'b','linewidth', 1);
legend('Conventional I-SGFD scheme','HEI-SGFD scheme','The difference between the two');
xlabel('Travel time/time step');
ylabel('Amp');
% ------------------------------
load('figure5_conventional_26ms.mat')
figure;plot(-Seismic_Txx(:,250),'r','linewidth', 1);
temp=Seismic_Txx;

load('figure5_HEI_26ms.mat')
hold on;plot(-Seismic_Txx(:,250),'k','linewidth', 1);
hold on;plot((-Seismic_Txx(:,250)+temp(:,250)),'b','linewidth', 1);
axis([0 2000 -15 10])
grid on
legend('Conventional I-SGFD scheme','HEI-SGFD scheme','The difference between the two');
xlabel('Travel time/time step');
ylabel('Amp');