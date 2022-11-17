clear;
clc
% close all

global v ratio M h tau

h=15;
tau=0.0035;
dt=tau;




r=v*dt/h;


tic;
coeffJune29=zeros(3600-1486+1,8);

options = optimset('Algorithm','trust-region-reflective','TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',4000,'MaxIter',4000);


i=1;
M=2;
x0=0.1*ones(1,M+2);
for v=(3600+50):-1:(3301+50)
    ratio=0.96;
    
    [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);   % Invoke optimizer
    x0=real(x);
    %        [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);
    %     x(end)=x(end)-abs(x(end))/10;
    coeffJune29(i,1:2)=x(1:2);
    coeffJune29(i,4:5)=x(3:4);
    %      coeffJune29(i,1:length(x))=x;
    
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

M=2;
x0=0.1*ones(1,M+2);
for v=(3300+50):-1:(3101+50)
    ratio=0.92;
    [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);   % Invoke optimizer
    %          [x,RESNORM,RESIDUAL,EXITFLAG]= lsqnonlin(@myfun7,x0);   % Invoke optimizer
    coeffJune29(i,1:2)=x(1:2);
    coeffJune29(i,4:5)=x(3:4);
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

M=3;
x0=0.1*ones(1,M+2);
for v=3100:-1:2801
    ratio=0.92;
    [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);   % Invoke optimizer
    %          [x,RESNORM,RESIDUAL,EXITFLAG]= lsqnonlin(@myfun7,x0);   % Invoke optimizer
    %     coeffJune29(i,1:length(x))=x;
    
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


for v=2800:-1:1486
    ratio=0.82;
    [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);   % Invoke optimizer
    %          [x,RESNORM,RESIDUAL,EXITFLAG]= lsqnonlin(@myfun7,x0);   % Invoke optimizer
    %     coeffJune29(i,1:length(x))=x;
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
save ('Conventional_35ms_Coeff_2.mat', 'coeffJune29');
