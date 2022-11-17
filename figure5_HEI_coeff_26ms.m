clear;
clc
% close all

global v ratio M h tau

h=15;
tau=0.0026;
M=3;

ratio=0.8;

x0=0.01*ones(1,M+2);%系数的初值是0

coeffJune=zeros(2300-857+1,5);
i=1;
tic
for v=1486/sqrt(3):3000/sqrt(3)
    
    ratio=0.8;
    [x,resnorm] = lsqnonlin(@myfun2,x0);   % Invoke optimizer
    coeffJune(i,:)=x;
    i=i+1;
end



ratio=0.5;
M=2;
ii=1;
x0=0.01*ones(1,M+2);%系数的初值是0
coeffJune2=zeros(800,5);
% for v=2300*sqrt(3):4790
for v=3000+600:4790+600
    [x,resnorm] = lsqnonlin(@myfun2,x0);   % Invoke optimizer
    coeffJune2(ii,1:2)=x(1:2);
    coeffJune2(ii,4:5)=x(3:4);
    ii=ii+1;
end

toc
% save('coeffJune.mat','coeffJune');
save('figure5_HEI_coeff_26ms.mat','coeffJune','coeffJune2');

