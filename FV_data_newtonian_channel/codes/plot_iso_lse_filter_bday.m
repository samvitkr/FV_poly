close all
clear


x1=100;
y1=100;
width=1100;
height=400;

nx=512;
nz=384;
ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
jcond=194;
%lt=0.05;
%lt=0.22;vim ca
%lt=-0.001;
lt=-0.005;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
mm=matfile('../data/mean_profiles.mat')
dUdy=reshape(mm.dUdy,[1,1,ny]);
dUdy=dUdy(1,1,111:end);
yp=yCheb(111:end)'+1;
[X,Z,Y]=meshgrid(xp,zp,yp);

f1=sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond');
m1=matfile(f1,'Writable',true);

ft =sprintf('../data/velgradfield_dfil_lseQ2_j_%03d.mat',jcond);
%ft =sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond);
%ft=sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond')
m=matfile(ft,'Writable',true)

ftu=sprintf('../data/velgradfield_ufil_lseQ2_j_%03d.mat',jcond);
%ftu=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond);
%ftu =sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond);
%ftu=sprintf('../data/velgrad_voz_field_lseQ4_j_156.mat')
mu=matfile(ftu,'Writable',true)

qrmsd=rms(m.lambda2,[1 2]);
qrmsu=rms(mu.lambda2,[1 2]);

% nld=(m.v).*(m.dvdx-m.dudy)-(m.w).*(m.dudz-m.dwdx);
% nlu=(mu.v).*(mu.dvdx-mu.dudy)-(mu.w).*(mu.dudz-mu.dwdx);
% syzd=nld+m.fx;
% syzu=nlu+mu.fx;
%nld=nld(:,:,111:end);
%nlu=nlu(:,:,111:end);
% % x1=150;
% % y1=150;
% % x2=1000;
% % y2=350;
 h1=figure('OuterPosition',...
    [x1 y1 width height]);
%%
subplot(1,3,1)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute(q,[2 1 3]), lt, permute((1e+3)*m.fxQ2,[2 1 3]))

isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(m1.lambda2,[2 1 3]), lt)%, permute(m.v,[2 1 3]))
colormap jet
%colorbar
lightangle(-60,-90)
axis equal
ylim([-0.8 1.2])
xlim([-0.5 0.5])
zlim([0 0.401])
yticks([-0.8:0.4:1.2])
zticks([0:0.2:0.4])

view(30,30)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
set(gca,'FontSize',12)
title('(a)')
%%
subplot(1,3,2)
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(m.lambda2,[2 1 3]), lt)%, permute(mu.v,[2 1 3])) 
colormap jet
%colorbar
lightangle(-60,-90)
axis equal
ylim([-0.8 1.2])
xlim([-0.5 0.5])
zlim([0 0.401])
yticks([-0.8:0.4:1.2])
zticks([0:0.2:0.4])
view(30,30)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
set(gca,'FontSize',12)
title('(b)')
%%

subplot(1,3,3)
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(mu.lambda2,[2 1 3]), lt)%, permute(mu.v,[2 1 3])) 
colormap jet
%colorbar
lightangle(-60,-90)
axis equal
ylim([-0.8 1.2])
xlim([-0.5 0.5])
zlim([0 0.401])
yticks([-0.8:0.4:1.2])
zticks([0:0.2:0.4])

view(30,30)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
set(gca,'FontSize',12)
title('(c)')

colormap parula
 %f1=sprintf("iso_l2tot_q2q4_%03d.fig",jcond)
%fn=sprintf("iso_l2Q2_fil_j_%03d.fig",jcond)
%saveas(h1,fn)
