clear;
clc
% close all

global v ratio M h tau

h=15;
tau=0.0035/2;
dt=tau;




r=v*dt/h;


tic;
coeffJune29=zeros(4790-1486+1,8);

options = optimset('Algorithm','trust-region-reflective','TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',4000,'MaxIter',4000);

i=1;
M=3;
x0=0.01*ones(1,M+2);

for v=1486:4790
    ratio=0.82;
    [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);   % Invoke optimizer
    coeffJune29(i,1:3)=x(1:3);
    coeffJune29(i,4:5)=x(4:5);
    x0=real(x);
    
    x2=real( x);
    temp=0;
    for ii=1:M
        temp=temp+2*x2(ii)*(-1)^(ii-1);
    end
    
    temp=temp-4*x2(M+1); %%%%%%%%%%%%%%%%
    temp=(1-4*x2(M+2))/temp;
    s2=sqrt(temp);
    coeffJune29(i,end-2)=s2;
    coeffJune29(i,end-1)=v*tau/h;
    %         coeffJune29(i,end-1)=v*tau/h;
    coeffJune29(i,end)= EXITFLAG;
    i=i+1;
    
end



toc
save ('Conventional_175ms_Coeff.mat', 'coeffJune29');
