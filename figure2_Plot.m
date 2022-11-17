
clear;clc
close all;
load('figure2_4ms_conventional.mat')

plotimage(-ptemp)
xlabel('x/dx')
ylabel('z/dz')
title('')
% dt
% h
% max(max(v))
% figure;plot(src)

load('figure2_4ms_HEI.mat')

plotimage(-ptemp)
xlabel('x/dx')
ylabel('z/dz')
% title('')
% dt
% h
% max(max(v))
% hold on;plot(src)

% plotimage(-ptemp3)
% xlabel('x/dx')
% ylabel('z/dz')
% title('')

% plotimage(-ptemp4)
% xlabel('x/dx')
% ylabel('z/dz')
% title('')

% plotimage(-p)
% xlabel('x/dx')
% ylabel('z/dz')
% title('')