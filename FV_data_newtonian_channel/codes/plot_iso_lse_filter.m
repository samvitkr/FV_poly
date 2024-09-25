close all
clear
nx=512;
nz=384;
ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
jcond=156;
%lt=0.05;
%lt=0.22;
lt=0.01;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
yp=yCheb(111:end)'+1;
[X,Z,Y]=meshgrid(xp,zp,yp);

ft =sprintf('../data/velgradfield_dfil_lseQ4_j_%03d.mat',jcond);
%ft =sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond);
%ft=sprintf('../data/velgrad_voz_field_lseQ2_j_156.mat')
m=matfile(ft,'Writable',true)

ftu=sprintf('../data/velgradfield_ufil_lseQ4_j_%03d.mat',jcond);
%ftu =sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond);
%ftu=sprintf('../data/velgrad_voz_field_lseQ4_j_156.mat')
mu=matfile(ftu,'Writable',true)

qrmsd=rms(m.Q,[1 2]);
qrmsu=rms(mu.Q,[1 2]);

nld=(m.v).*(m.dvdx-m.dudz)-(m.w).*(m.dudz-m.dwdx);
nlu=(mu.v).*(mu.dvdx-mu.dudz)-(mu.w).*(mu.dudz-mu.dwdx);
syzd=nld+m.fx;
syzu=nlu+mu.fx;
%nld=nld(:,:,111:end);
%nlu=nlu(:,:,111:end);
x1=150;
y1=150;
x2=1000;
y2=350;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
subplot(1,2,1)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute(q,[2 1 3]), lt, permute((1e+3)*m.fxQ2,[2 1 3]))

isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(m.Q,[2 1 3]), lt, permute(-syzd,[2 1 3]))

%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fxQ2,[2 1 3]), -fth)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fxQ2,[2 1 3]), fth)
%clim([-1 1])
colormap jet
colorbar
lightangle(-45,-90)
axis equal
% ylim([-1 1])
% xlim([-0.5 0.5])
% zlim([0 0.4])
view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on

subplot(1,2,2)
%isosurface( permute(mu.Z,[2 1 3]), permute(mu.X,[2 1 3]), permute(mu.Y,[2 1 3]),...
%permute(qu,[2 1 3]), lt, permute((1e+3)*mu.fxQ2,[2 1 3])) 

isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(mu.Q,[2 1 3]), lt, permute(-syzu,[2 1 3])) 

%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fyQ2,[2 1 3]), -fth)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fyQ2,[2 1 3]), fth)
%clim([-1 1])
colormap jet
colorbar
lightangle(-45,-90)
axis equal
% ylim([-1 1])
% xlim([-0.5 0.5])
% zlim([0 0.4])
view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on


% f1=sprintf("iso_lse_vel_nl_Q5e-3_filter_j2_%03d.fig",jcond)
% saveas(h1,f1)
% fl=sprintf("lse_profiles_j_%03d.mat",jcond)
% mfl=matfile(fl,'Writable',true)
% mfl.lt=lt;
% 
% mfl.voz2=squeeze(mean(voz2,[1 2]));
% mfl.voz4=squeeze(mean(voz4,[1 2]));
% mfl.woy2=squeeze(mean(woy2,[1 2]));
% mfl.woy4=squeeze(mean(woy4,[1 2]));
% 
% mfl.voz2u=squeeze(mean(voz2u,[1 2]));
% mfl.voz4u=squeeze(mean(voz4u,[1 2]));
% mfl.woy2u=squeeze(mean(woy2u,[1 2]));
% mfl.woy4u=squeeze(mean(woy4u,[1 2]));
% 
% mfl.voz2q=squeeze(mean(voz2.*(q2>lt),[1 2]));
% mfl.voz4q=squeeze(mean(voz4.*(q4>lt),[1 2]));
% mfl.woy2q=squeeze(mean(woy2.*(q2>lt),[1 2]));
% mfl.woy4q=squeeze(mean(woy4.*(q4>lt),[1 2]));
% 
% mfl.voz2uq=squeeze(mean(voz2u.*(qu2>lt),[1 2]));
% mfl.voz4uq=squeeze(mean(voz4u.*(qu4>lt),[1 2]));
% mfl.woy2uq=squeeze(mean(woy2u.*(qu2>lt),[1 2]));
% mfl.woy4uq=squeeze(mean(woy4u.*(qu4>lt),[1 2]));
% 
% mfl.fx2=squeeze(mean(m.fxQ2,[1 2]));
% mfl.fx4=squeeze(mean(m.fxQ4,[1 2]));
% mfl.fx2u=squeeze(mean(mu.fxQ2,[1 2]));
% mfl.fx4u=squeeze(mean(mu.fxQ4,[1 2]));
% 
% mfl.fx2q=squeeze(mean(m.fxQ2.*(q2>lt),[1 2]));
% mfl.fx4q=squeeze(mean(m.fxQ4.*(q4>lt),[1 2]));
% mfl.fx2uq=squeeze(mean(mu.fxQ2.*(qu2>lt),[1 2]));
% mfl.fx4uq=squeeze(mean(mu.fxQ4.*(qu4>lt),[1 2]));


