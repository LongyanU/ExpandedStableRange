clear;
clc
% close all

global v ratio M h tau

h=15;
tau=0.0035;
dt=tau;




r=v*dt/h;


tic;
coeffJune29=zeros(4790-1486+1,8);

options = optimset('Algorithm','trust-region-reflective','TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',4000,'MaxIter',4000);

i=1;
M=3;
x0=0.01*ones(1,M+2);

for v=1486:3000
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


M=2;
x0=0.1*ones(1,M+2);
for v=3001:3500
    ratio=0.96;
    [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);   % Invoke optimizer
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
    coeffJune29(i,end)= EXITFLAG;
    i=i+1;
    
end

M=2;
% x0=0.1*ones(1,M+2);
for v=3501+50:3600+50
    ratio=0.96;
    [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);   % Invoke optimizer

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
    coeffJune29(i,end)= EXITFLAG;
    i=i+1;
    
end

M=3;
x0=0.1*ones(1,M+2);
for v=3601+200:4200+200
    ratio=0.999;
    [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);   % Invoke optimizer

%     coeffJune29(i,1:2)=x(1:2);
%     coeffJune29(i,4:5)=x(3:4);
    
     coeffJune29(i,1:5)=x(1:5);
   
    
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
    coeffJune29(i,end)= EXITFLAG;
    i=i+1;
    
end

M=3;
x0=0.1*ones(1,M+2);
for v=4201+250:4400+250
    ratio=0.9;
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


M=3;
x0=0.01*ones(1,M+2);
for v=4401+300:4500+300
    ratio=0.9;
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
% 
M=3;
x0=0.01*ones(1,M+2);
for v=4501+400:4600+400
    ratio=0.92;
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

M=2;
x0=0.01*ones(1,M+2);
for v=4601+550:4800+550
    ratio=0.92;
    [x,RESNORM,RESIDUAL,EXITFLAG] = lsqnonlin(@myfun7,x0);   % Invoke optimizer

%     coeffJune29(i,1:3)=x(1:3);
%     coeffJune29(i,4:5)=x(4:5);
    
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

    coeffJune29(i,end)= EXITFLAG;
    i=i+1;
    
end

toc
save ('Conventional_35ms_Coeff.mat', 'coeffJune29');
