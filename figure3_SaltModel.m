clear;clc;close all



load('figure3_175ms_Conventional.mat')
figure;imagesc(v(45:end-45,90:end-90))
xlabel('x/dx')
ylabel('z/dz')

h=colorbar;
set(get(h,'title'),'string','m/s');
set(gca,'FontSize',18,'FontName','Times New Roman')

% subplot(2,2,1)
% imagesc(seis_recordp(1:1250*2,80:nx-80),[-3.5*10 7.3*10]);
% colormap gray
% % xlabel('x/dx')
% ylabel('Travel time/time step')
% text(-100,-40,'(a)')
% 
% 
% load('figure3_175ms_HEI.mat')
% subplot(2,2,2)
% imagesc(seis_recordp(1:1250*2,80:nx-80),[-3.5*10 7.3*10]);
% colormap gray
% % xlabel('x/dx')
% % ylabel('Travel time/time step')
% text(-100,-40,'(b)')
% 
% load('figure3_35ms_Conventional.mat')
% subplot(2,2,3)
% imagesc(seis_recordp(1:1250,80:nx-80),[-3.5*10 7.3*10]);
% colormap gray
% xlabel('x/dx')
% ylabel('Travel time/time step')
% text(-100,-40,'(c)')
% 
% 
% load('figure3_35ms_HEI.mat')
% subplot(2,2,4)
% imagesc(seis_recordp(1:1250,80:nx-80),[-3.5*10 7.3*10]);
% colormap gray
% xlabel('x/dx')
% % ylabel('Travel time/time step')
% text(-100,-40,'(d)')
% 
% 
