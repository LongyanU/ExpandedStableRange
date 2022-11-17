clear;clc;close all
load('figure3_175ms_Conventional.mat')
plotimage(-seis_recordp(1:1250*2,80:nx-80))
xlabel('x/dx')
ylabel('Travel time/time step')
title('')
% % tempaa=seis_recordp;
% % figure;plot(-seis_recordp(1:2500,380)./max(-seis_recordp(1:2500,380)),'r','linewidth',1);


load('figure3_175ms_HEI.mat')
plotimage(-seis_recordp(1:1250*2,80:nx-80))
xlabel('x/dx')
ylabel('Travel time/time step')
title('')

% % hold on;plot(-seis_recordp(1:2500,380)./max(-seis_recordp(1:2500,380)),'k','linewidth',1);
% % hold on;plot(seis_recordp(1:2500,380)./max(-seis_recordp(1:2500,380))-tempaa(1:2500,380)./max(-tempaa(1:2500,380)),'b','linewidth',1);
% % axis([1 2500 -0.53 1.001])
% % grid on
% % legend('Conventional I-SGFD scheme','HEI-SGFD scheme','The difference between the two');
% % xlabel('Travel time/time step')
% % ylabel('Amp')