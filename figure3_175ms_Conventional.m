% 历时 1215.737671 秒。2520*2
% 历时 604.071077 秒。
clear
clc %%%%%%%
close all
nt=1250*2;   % number of time steps
eps=.6;     % stability
isnap=10;    % snapshot sampling


load('vv')

c1=flipud(c);

v=c1;
nx=799;
nx=nx+45*4;
nz=475;
nz=nz+45*2;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90*2
        vv(ii+45,jj+45)=v(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45*2
        vv(ii+45,jj)=v(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45*4:nx
        vv(ii+45,jj)=v(ii,800);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end


% % % % % clear v
% % % % %
% % % % %
% % % % % v=vv;
% % % % % [nz, nx]=size(v);

clear v
v=vv;

% % % % itemp=1;
% % % % jtemp=1;
% % % % for i=2:2:nz
% % % %     for j=2:2:nx
% % % %         vtemp(floor(i/2),floor(j/2))=v(i,j);
% % % %     end
% % % %     %      jtemp=jtemp+1;
% % % %     %     itemp=itemp+1;
% % % % end
% % % % 
% % % % clear v
% % % % v=vtemp;
% % % % [nz,nx]=size(vtemp);


%
% % % % for ii=1:nz
% % % %     for jj=1:nx
% % % %         %         if(v(ii,jj)) <=2250
% % % %         %             v(ii,jj)=2250;
% % % %         %         end
% % % %         
% % % % %         if(v(ii,jj)) >3600
% % % % %             v(ii,jj)=3600;
% % % % %         end
% % % %         
% % % %     end
% % % % end


dx=15;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.0035/2; % calculate time step from stability criterion
tau=dt;


f0=40;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^8*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=(diff(src))/dx^2;				% time derivative to obtain gaussian


zs=46;
xs=600-150+25;
% xs=floor((600-150+25)/2);
xs=floor(nx/2);

seis_recordVx=zeros(nt,nx);
seis_recordVz=zeros(nt,nx);
seis_recordp=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;

VxzAddedpoint=p;

r=v*dt/h;

M=3;
load ('Conventional_175ms_Coeff.mat')
coeffJune29=real(coeffJune29);

coeff=zeros(nz,nx,5);
for ii=1:nz
    for jj=1:nx
%         if ( v(ii,jj)<=3600 )
            
            coeff(ii,jj,:)=real(coeffJune29((floor(abs(v(ii,jj)-1486)))+1,1:5));
            
%         else
            %
%             coeff(ii,jj,:)= [0.5675    0.1101   -0.0027    0.0569    0.0877];
%                          coeff(ii,jj,:)= [0.504078, 0.146255,0, 0.0272974, 0.156212];
            % stable however not accurate
            % v=5000;
            %
            % h=15;
            % tau=0.0035;
            % M=3;
            % r=v*tau/h;
            % ratio=0.8;
            % stable however not accurate
%         end
    end
end

b=1-2*coeff(:,:,end);
a = coeff(:,:,end);
cc = a;

for j=1:nx,
    A1{j} = (gallery('tridiag',a(1:nz-1,j),b(:,j),cc(1:nz-1,j)));
end



A2=cell(nz,1);
for k=1:nz,
    A2 {k}= (gallery('tridiag',a(k,1:nx-1),b(k,:),cc(k,1:nx-1) ));      %%解三对角阵，直接matlab解了
end


Vx=zeros(nz,nx);
Vz=zeros(nz,nx);
d2pzz=p;d2pxx=p;


tic
for it=1:nt-2,
    
    
    d2px11=(circshift(Vx,[0 0])-circshift(Vx,[0 1]));
    d2px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
    d2px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
    
    addPoint=(circshift(Vx,[-1 0])-circshift(Vx,[-1 1]));
    addPoint=addPoint+(circshift(Vx,[1 0])-circshift(Vx,[1 1]));
    d2px=coeff(:,:,1).*d2px11+coeff(:,:,2).*d2px12+coeff(:,:,3).*d2px13;
    d2px=d2px+addPoint.*coeff(:,:,4);
    
    d2pz11=(circshift(Vz,[0 0])-circshift(Vz,[-1 0]));
    d2pz12=(circshift(Vz,[1 0])-circshift(Vz,[-2 0]));
    d2pz13=(circshift(Vz,[2 0])-circshift(Vz,[-3 0]));
    
    addPoint=(circshift(Vz,[0 -1])-circshift(Vz,[-1 -1]));
    addPoint=addPoint+(circshift(Vz,[0 1]) -circshift(Vz,[-1 1]));
    
    d2pz=coeff(:,:,1).*d2pz11+coeff(:,:,2).*d2pz12+coeff(:,:,3).*d2pz13;
    d2pz=d2pz+addPoint.*coeff(:,:,4);
    
    for j=1:nx,
        d2pzz(:,j)=A1{j}\d2pz(:,j);
    end
    for k=1:nz,
        d2pxx(k,:)=A2{k}\d2px(k,:)';
    end
    p=p-dt*v.^2.*(d2pxx+d2pzz)/h;
    p(zs,xs)=p(zs,xs)+(src(it)+src(it+1)/2);
    %     [p,p]=spongeABC(p,p,nx,nz,45,45,0.009);
    % time saved
    seis_recordp(it,:)=p(zs,:);
    
    d2px1=(circshift(p,[0 -1])-circshift(p,[0 0]));
    d2px2=(circshift(p,[0 -2])-circshift(p,[0 1]));
    d2px3=(circshift(p,[0 -3])-circshift(p,[0 2]));
    addPoint=(circshift(p,[-1 -1])-circshift(p,[-1 0]));
    addPoint=addPoint+(circshift(p,[1 -1])-circshift(p,[1 0]));
    d2px=coeff(:,:,1).*d2px1+coeff(:,:,2).*d2px2+coeff(:,:,3).*d2px3;
    d2px=d2px+addPoint.*coeff(:,:,4);
    
    d2pz1=(circshift(p,[1])-circshift(p,[0]));
    d2pz2=(circshift(p,[2])-circshift(p,[-1]));
    d2pz3=(circshift(p,[3])-circshift(p,[-2]));
    addPoint=(circshift(p,[1 -1])-circshift(p,[0 -1]));
    addPoint=addPoint+(circshift(p,[1 1]) -circshift(p,[0 1]));
    d2pz=coeff(:,:,1).*d2pz1+coeff(:,:,2).*d2pz2+coeff(:,:,3).*d2pz3;
    d2pz=d2pz+addPoint.*coeff(:,:,4);
    
    for j=1:nx,
        d2pz1(:,j)=A1{j}\d2pz(:,j);
    end
    for k=1:nz,
        d2px1(k,:)=A2{k}\d2px(k,:)';
    end
    
    
    Vx=Vx-dt*d2px1/h;
    Vz=Vz-dt*d2pz1/h;
    
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,85,45,0.009);
    
    seis_recordVx(it,:)=Vx(zs,:);
    seis_recordVz(it,:)=Vz(zs,:);
    
    if rem(it,isnap)== 0,
        imagesc(x,z,real(p)), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(real(p)))))
        drawnow
    end
end
toc
save('figure3_175ms_Conventional.mat')
