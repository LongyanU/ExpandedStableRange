clear;
clc
% close all

global v ratio M h tau

h=15;
tau=0.0026/2;
M=3;

ratio=0.8;

x0=0.01*ones(1,M+2);%系数的初值是0

coeffJune=zeros(ceil(4790/sqrt(3))-857+1,5);
i=1;
tic
for v=1486/sqrt(3):4790/sqrt(3)
    
    ratio=0.8;
    [x,resnorm] = lsqnonlin(@myfun7,x0);   % Invoke optimizer
    coeffJune(i,:)=x;
    i=i+1;
end

toc
% save('coeffJune.mat','coeffJune');
save('figure5_conventional_13ms_coeff.mat','coeffJune');

