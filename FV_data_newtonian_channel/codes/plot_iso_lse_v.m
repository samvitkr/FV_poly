close all
clear
nx=512;
nz=384;
ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
jcond=145;
%lt=0.05;
%lt=0.22;vim ca
lt=-0.005;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
mm=matfile('../data/mean_profiles.mat')
ut=mm.ret/mm.re;
dUdy=reshape(mm.dUdy,[1,1,ny]);
dUdy=dUdy(1,1,111:end);
yp=yCheb(111:end)'+1;
[X,Z,Y]=meshgrid(xp,zp,yp);

%ft =sprintf('../data/velgradfield_dfil_lseQ4_j_%03d.mat',jcond);
ft =sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond);
%ft=sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond')
m=matfile(ft,'Writable',true)
ftu =sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond);
%ftu=sprintf('../data/velgrad_voz_field_lseQ4_j_156.mat')
%ftu=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond')

mu=matfile(ftu,'Writable',true)


vozd=(m.v).*(m.dvdx-m.dudy-dUdy);
woyd=(m.w).*(m.dudz-m.dwdx);
vozu=(mu.v).*(mu.dvdx-mu.dudy-dUdy);
woyu=(mu.w).*(mu.dudz-mu.dwdx);

x1=150;
y1=150;
x2=400;
y2=400;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
subplot(2,2,1)
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(m.lambda2,[2 1 3]), lt, permute(vozd./(-ut^2),[2 1 3]))

c1=colorbar;
clim([-1 1])
ylabel(c1,"$ v\omega_z/(-u_{\tau^2}/H) $",'interpreter','latex')

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

subplot(2,2,2)
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(mu.lambda2,[2 1 3]), lt, permute(vozu./(-ut^2),[2 1 3])) 

c2=colorbar;
clim([-1 1])

ylabel(c2,"$ v\omega_z/(-u_{\tau^2}/H) $",'interpreter','latex')
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

subplot(2,2,3)
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(m.lambda2,[2 1 3]), lt, permute(-woyd./(-ut^2),[2 1 3]))

c3=colorbar;
clim([-1 1])

ylabel(c3,"$ -w\omega_y/(-u_{\tau^2}/H) $",'interpreter','latex')
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

subplot(2,2,4)
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(mu.lambda2,[2 1 3]), lt, permute(-woyu./(-ut^2),[2 1 3])) 
colormap redblue
c4=colorbar;
clim([-1 1])

ylabel(c4,"$ -w\omega_y/(-u_{\tau^2}/H) $",'interpreter','latex')
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

%f1=sprintf("iso_l2_vpn_nl_%03d.fig",jcond)
 %f1=sprintf("iso_l2_vpn_%03d.fig",jcond)
%saveas(h1,f1)
