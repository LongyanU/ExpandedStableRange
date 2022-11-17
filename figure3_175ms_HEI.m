
% 历时 312.524848 秒。
clear
clc %%%%%%%
close all
nt=1250*2;    % number of time steps
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


clear v


v=vv;
[nz, nx]=size(v);

% % % % % clear v
% % % % % v=vv;

% % % % % % itemp=1;
% % % % % % jtemp=1;
% % % % % for i=2:2:nz
% % % % %     for j=2:2:nx
% % % % %         vtemp(floor(i/2),floor(j/2))=v(i,j);
% % % % %         
% % % % %     end
% % % % %     %      jtemp=jtemp+1;
% % % % %     %     itemp=itemp+1;
% % % % % end
% % % % % 
% % % % % clear v
% % % % % v=vtemp;
% % % % % [nz,nx]=size(vtemp);


%
% % % % % % for ii=1:nz
% % % % % %     for jj=1:nx
% % % % % %         if(v(ii,jj)) >=4700
% % % % % %             v(ii,jj)=4700;
% % % % % %         end
% % % % % %         
% % % % % %         %         if(v(ii,jj)) >4148
% % % % % %         %             v(ii,jj)=4148;
% % % % % %         %         end
% % % % % %         
% % % % % %     end
% % % % % % end


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
load ('figure3_175ms_HEI_Coeff.mat')
coeffJune29=real(coeffJune29);

coeff=zeros(nz,nx,5);
for ii=1:nz
    for jj=1:nx
%         if ( v(ii,jj)>4600 && v(ii,jj)<4791 ) 
%             coeff(ii,jj,:)=[0.840434, -0.0369836, 0, 0.188318, -0.32232];
%         elseif ( v(ii,jj)>4500 &&v(ii,jj)<=4600 )  
%             coeff(ii,jj,:)=[0.590719, 0.0564095,0, 0.142804, -0.144154];
%         elseif ( v(ii,jj)>4400 && v(ii,jj)<=4500 )  
%             coeff(ii,jj,:)=[0.58087, 0.0644378,0, 0.135598, -0.113558];
%         elseif ( v(ii,jj)>4300 && v(ii,jj)<=4400 ) 
%             coeff(ii,jj,:)=[0.698613, 0.0335653, -0.034505, 0.188905, -0.184008];
%         else
            coeff(ii,jj,:)=real(coeffJune29(abs(floor((v(ii,jj)-1486)))+1,1:5));
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
    
    addedShifitedMatrix=zeros(nz,nx);
    for ii=1:M
        addedShifitedMatrix=addedShifitedMatrix+ coeff(:,:,ii).*(circshift(Vx,[0 1-ii])-circshift(Vx,[0 ii])) ;
    end
    
    addPoint=(circshift(Vx,[-1 0])-circshift(Vx,[-1 1]));
    addPoint=addPoint+(circshift(Vx,[1 0])-circshift(Vx,[1 1]));
    d2px=addedShifitedMatrix+addPoint.*coeff(:,:,M+1);
    
    addedShifitedMatrix=zeros(nz,nx);
    for ii=1:M
        addedShifitedMatrix=addedShifitedMatrix+ coeff(:,:,ii).*(circshift(Vz,[ii-1 0])-circshift(Vz,[-ii 0])) ;
    end
    
    
    addPoint=(circshift(Vz,[0 -1])-circshift(Vz,[-1 -1]));
    addPoint=addPoint+(circshift(Vz,[0 1]) -circshift(Vz,[-1 1]));
    
    d2pz=addedShifitedMatrix+addPoint.*coeff(:,:,M+1);
    
    for j=1:nx,
        d2pzz(:,j)=A1{j}\d2pz(:,j);
    end
    for k=1:nz,
        d2pxx(k,:)=A2{k}\d2px(k,:)';
    end
    p=p-dt*v.^2.*(d2pxx+d2pzz)/h;
    p(zs,xs)=p(zs,xs)+src(it);
    seis_recordp(it,:)=p(zs,:);
    
    % time saved
    d2px1=circshift(p,[0 -1])-circshift(p,[0 0]);
    d2pz1=circshift(p,[1])-p;
    
    Vx=Vx-dt*d2px1/h;
    Vz=Vz-dt*d2pz1/h;
    
%     [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,85,45,0.007);
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
save('figure3_175ms_HEI.mat')
