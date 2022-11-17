clear
clc %%%%%%%
close all
nt=2002;    % number of time steps
eps=.6;     % stability
isnap=20;    % snapshot sampling
load('vv')
% 历时 427.617723 秒。 Nov 12 ,2022
c1=flipud(c);

v=c1;
nx=800;
nx=nx+45*2;
nz=475;
nz=nz+45*2;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=v(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=v(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
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




clear v
v=vv;

% for i=2:2:nz
%     for j=2:2:nx
%         vtemp(floor(i/2),floor(j/2))=v(i,j);
%     end
% end
% 
% clear v
% v=vtemp;
[nz,nx]=size(v);
% 
% for i=1:nz
%     for j=1:nx
%         if(v(i,j)>4000)
%             v(i,j)=4000;
%         end
%     end
% end

vp=v;
vs=vp/sqrt(3);
rho=1000*ones(nz,nx);

dx=15;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.0026; % calculate time step from stability criterion
tau=dt;


f0=35;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian


xs=floor(nx/2);
zs=46;

seis_recordVx=zeros(nt,nx);
seis_recordVz=zeros(nt,nx);
seis_recordTxx=zeros(nt,nx);
seis_recordTzz=zeros(nt,nx);
seis_recordTxz=zeros(nt,nx);

%  ******************************************************************
%  txz-------------vx
%  |                |
%  |                |
%  |                |
%  |                |
%  |                |
%  vz---------------txx,tzz
% Virieux 1986 Geophysics

dx=h;
dz=h;

p=zeros([nz nx]); Vx=p; Vz=p;
Txxx=p;
Txzz=p;
Tzzz=p;
Txzx=p;

Txx=p;
Txz=p;
Tzz=p;


Vxx=p;
Vzz=p;
Vxz=p;
Vzx=p;

temp1=min(min(vs));
load('figure5_HEI_coeff_26ms.mat')
coeffJune=real(coeffJune);
coeffJune2=real(coeffJune2);
coeff=zeros(nz,nx,5);
for ii=1:nz
    for jj=1:nx
        
        if (vp(ii,jj)<=3000)
            coeff(ii,jj,:)=real(coeffJune(floor(vs(ii,jj)-temp1)+1,:));
        else
            %             coeff(ii,jj,:)=coeffJune2(floor(vp(ii,jj)-2300*sqrt(3))+1,:);
            coeff(ii,jj,:)=coeffJune2(floor(vp(ii,jj)-3000)+1,:);
        end
        
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




Seismic_Vx=zeros(nt,nx);
Seismic_Vz=zeros(nt,nx);
Seismic_Txx=zeros(nt,nx);
Seismic_Tzz=zeros(nt,nx);
Seismic_Txz=zeros(nt,nx);

tic
for it=1:nt-3,
    
    Txxx=coeff(:,:,1).*( (Txx)-circshift(Txx,[0,-1]))+...
        coeff(:,:,2).*(circshift(Txx,[0,1])-circshift(Txx,[0,-2]))+...
        coeff(:,:,3).*( circshift(Txx,[0,2])-circshift(Txx,[0,-3]));
    
    addPoint=(circshift(Txx,[-1 0])-circshift(Txx,[-1 -1]));
    addPoint=addPoint+(circshift(Txx,[1 0])-circshift(Txx,[1 -1]));
    
    Txxx=addPoint.*coeff(:,:,4)+Txxx;
    
    for k=1:nz,
        Txxx(k,:)=A2{k}\Txxx(k,:)';
    end
    
    %Txz/z
    Txzz=circshift(Txz,[ 1])-circshift(Txz,[ 0]);
    %     Txzz=coeff(:,:,1).*(circshift(Txz,[ 1])-circshift(Txz,[ 0]))+...
    %         coeff(:,:,2).*(circshift(Txz,[ 2])-circshift(Txz,[ -1]))+...
    %         coeff(:,:,3).*(circshift(Txz,[ 3])-circshift(Txz,[ -2]));
    %
    %     addPoint=circshift(Txz,[ 1  -1])-circshift(Txz,[ 0  -1]);
    %     addPoint=addPoint+(circshift(Txz,[ 1  1])-circshift(Txz,[ 0  1]));
    %
    %     Txzz=addPoint.*coeff(:,:,4)+Txzz;
    %     for j=1:nx,
    %         Txzz(:,j)=A1{j}\Txzz(:,j);
    %     end
    
    
    %Tzz/z
    Tzzz=coeff(:,:,1).*( (Tzz)-circshift(Tzz,[-1]))+...
        coeff(:,:,2).*(circshift(Tzz,[1])-circshift(Tzz,[-2]))+...
        coeff(:,:,3).*( circshift(Tzz,[2])-circshift(Tzz,[-3]));
    
    addPoint=(circshift(Tzz,[ 0  -1])-circshift(Tzz,[ -1  -1]));
    addPoint=addPoint+(circshift(Tzz,[ 0  1])-circshift(Tzz,[ -1  1]));
    
    Tzzz=addPoint.*coeff(:,:,4)+Tzzz;
    for j=1:nx,
        Tzzz(:,j)=A1{j}\Tzzz(:,j);
    end
    
    %Txz/x
    Txzx=circshift(Txz,[0 1])-Txz;
    %     Txzx=coeff(:,:,1).*(circshift(Txz,[0 1])-Txz)+...
    %         coeff(:,:,2).*(circshift(Txz,[0 2])-circshift(Txz,[0 -1]))+...
    %         coeff(:,:,3).*(circshift(Txz,[0 3])-circshift(Txz,[0 -2]));
    %
    %     addPoint=(circshift(Txz,[-1 1])-circshift(Txz,[-1 0]));
    %     addPoint=addPoint+(circshift(Txz,[1 1])-circshift(Txz,[1 0]));
    %
    %     Txzx=addPoint.*coeff(:,:,4)+Txzx;
    %
    %     for k=1:nz,
    %         Txzx(k,:)=A2{k}\Txzx(k,:)';
    %     end
    
    Vx=Vx+1./(rho).*dt.*(Txxx+Txzz)/h;
    Vz=Vz+1./(rho).*dt.*(Tzzz+Txzx)/h;
    
    Seismic_Vx(it,:)=Vx(46,:);
    Seismic_Vz(it,:)=Vz(46,:);
    
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    Vxz=coeff(:,:,1).*(circshift(Vx,[0])-circshift(Vx,[-1]))+...
        coeff(:,:,2).*(circshift(Vx,[1])-circshift(Vx,[-2]))+...
        coeff(:,:,3).*(circshift(Vx,[2])-circshift(Vx,[-3]));
    
    addPoint=circshift(Vx,[ 0  -1])-circshift(Vx,[ -1  -1]);
    addPoint=addPoint+(circshift(Vx,[ 0  1])-circshift(Vx,[ -1  1]));
    
    Vxz=addPoint.*coeff(:,:,4)+Vxz;
    for j=1:nx,
        Vxz(:,j)=A1{j}\Vxz(:,j);
    end
    
    Vzz=(circshift(Vz,[1])-circshift(Vz,[0]));
    %     Vzz=coeff(:,:,1).*(circshift(Vz,[1])-circshift(Vz,[0]))+...
    %         coeff(:,:,2).*(circshift(Vz,[2])-circshift(Vz,[-1]))+...
    %         coeff(:,:,3).*(circshift(Vz,[3])-circshift(Vz,[-2]));
    %
    %     addPoint=circshift(Vz,[ 1  -1])-circshift(Vz,[ 0  -1]);
    %     addPoint=addPoint+(circshift(Vz,[ 1  1])-circshift(Vz,[ 0  1]));
    %
    %     Vzz=addPoint.*coeff(:,:,4)+Vzz;
    %     for j=1:nx,
    %         Vzz(:,j)=A1{j}\Vzz(:,j);
    %     end
    
    Vxx=(circshift(Vx,[0 1])-circshift(Vx,[0 0]));
    %     Vxx=coeff(:,:,1).*(circshift(Vx,[0 1])-circshift(Vx,[0 0]))+...
    %         coeff(:,:,2).*(circshift(Vx,[0 2])-circshift(Vx,[0 -1]))+...
    %         coeff(:,:,3).*(circshift(Vx,[0 3])-circshift(Vx,[0 -2]));
    %
    %     addPoint=(circshift(Vx,[-1 1])-circshift(Vx,[-1 0]));
    %     addPoint=addPoint+(circshift(Vx,[1 1])-circshift(Vx,[1 0]));
    %
    %     Vxx=addPoint.*coeff(:,:,4)+Vxx;
    %
    %     for k=1:nz,
    %         Vxx(k,:)=A2{k}\Vxx(k,:)';
    %     end
    %
    
    Vzx=coeff(:,:,1).*(circshift(Vz,[0 0])-circshift(Vz,[0 -1]))+...
        coeff(:,:,2).*(circshift(Vz,[0 1])-circshift(Vz,[0 -2]))+...
        coeff(:,:,3).*(circshift(Vz,[0 2])-circshift(Vz,[0 -3]));
    
    addPoint=(circshift(Vz,[-1 0])-circshift(Vz,[-1 -1]));
    addPoint=addPoint+(circshift(Vz,[1 0])-circshift(Vz,[1 -1]));
    
    Vzx=addPoint.*coeff(:,:,4)+Vzx;
    
    for k=1:nz,
        Vzx(k,:)=A2{k}\Vzx(k,:)';
    end
    
    Txx=Txx+dt*rho.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*rho.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*rho.*(vs.^2).*(Vxz+Vzx)/h;
    
    Txx(zs,xs)=Txx(zs,xs)+src(it);
    Tzz(zs,xs)=Tzz(zs,xs)+src(it);
    
    Seismic_Txx(it,:)=Txx(46,:);
    Seismic_Tzz(it,:)=Tzz(46,:);
    Seismic_Txz(it,:)=Txz(46,:);
    
    if rem(it,isnap)== 0,
        %         imagesc(-Vx(46*1:end,:),[-10^5 10^5]), axis equal
        imagesc(-Txx(46:end,:)), axis equal
        %         imagesc(-Vx,[-10^5 10^5]), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(Vx))))
        drawnow
    end
end
toc

save('figure5_HEI_26ms.mat')
