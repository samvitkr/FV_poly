close all
clear


x1=100;
y1=100;
width=1400;
height=450;

nx=512;
nz=384;
ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
jcond=194;
jc=jcond-110;
%lt=0.05;
%lt=0.22;vim ca
%lt=-0.001;
lt=-0.01;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
load('luldrms.mat');
mm=matfile('../data/mean_profiles.mat')
dUdy=reshape(mm.dUdy,[1,1,ny]);
dUdy=dUdy(1,1,111:end);
yp=yCheb(111:end)'+1;
[X,Z,Y]=meshgrid(xp,zp,yp);

%f1=sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond');
%m1=matfile(f1,'Writable',true);

%ft =sprintf('../data/velgradfield_dfil_lseQ2_j_%03d.mat',jcond);
%ft =sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond);
ft=sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond)
m2=matfile(ft,'Writable',true)

%ftu=sprintf('../data/velgradfield_ufil_lseQ4_j_%03d.mat',jcond);
ftu=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond);
%ftu =sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond);
%ftu=sprintf('../data/velgrad_voz_field_lseQ4_j_156.mat')
m4=matfile(ftu,'Writable',true)

%qrmsd=rms(m.lambda2,[1 2]);
%qrmsu=rms(mu.lambda2,[1 2]);

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
m2lrms=squeeze(rms(m2.lambda2,[1 2]));
m2lt=-m2lrms(jc);
m2.lt=m2lt;
subplot(1,2,1)
oz=m2.dvdx-m2.dudy;%-m2.dUdym;
oy=m2.dudz-m2.dwdx;
ox=m2.dwdy-m2.dvdz;
omag=sqrt(ox.^2+oy.^2+oz.^2);
m2cos=oz./omag;
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(m2.lambda2,[2 1 3]), m2lt, permute(m2cos,[2 1 3]))
% colormap redblue
clim([-1 1])
%colorbar
lightangle(-60,-90)
axis equal
% xlim([-0.301 0.301])
% ylim([-0.8 1.01])
% zlim([0 0.301])
% yticks([-0.8:0.2:1])
% zticks([0:0.3:0.3])
% xticks([-0.3:0.3:0.3])

% xlim([-0.301 0.301])
% ylim([-0.8 0.6])
% zlim([0 0.401])
% yticks([-0.8:0.2:0.6])
% zticks([0:0.2:0.4])
% xticks([-0.3:0.3:0.3])

% xlim([-0.501 0.501])
% ylim([-1.2 0.6])
% zlim([0 0.801])
% yticks([-1.2:0.3:0.6])
% zticks([0:0.4:0.8])
% xticks([-0.5:0.5:0.5])

view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
set(gca,'FontSize',12)
title('(a)')
%%

m4lrms=squeeze(rms(m4.lambda2,[1 2]));
m4lt=-m4lrms(jc);
m4.lt=m4lt;
oz=m4.dvdx-m4.dudy;%-m4.dUdym;
oy=m4.dudz-m4.dwdx;
ox=m4.dwdy-m4.dvdz;
omag=sqrt(ox.^2+oy.^2+oz.^2);
m4cos=oz./omag;
subplot(1,2,2)
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(m4.lambda2,[2 1 3]), m4lt, permute(m4cos,[2 1 3])) 
colormap redblue
clim([-1 1])
c=colorbar;
lightangle(-60,-90)
axis equal

% xlim([-0.301 0.301])
% ylim([-0.8 1.01])
% zlim([0 0.301])
% yticks([-0.8:0.2:1])
% zticks([0:0.3:0.3])
% xticks([-0.3:0.3:0.3])

% xlim([-0.301 0.301])
% ylim([-0.8 0.6])
% zlim([0 0.401])
% yticks([-0.8:0.2:0.6])
% zticks([0:0.2:0.4])
% xticks([-0.3:0.3:0.3])

% xlim([-0.501 0.501])
% ylim([-1.2 0.6])
% zlim([0 0.801])
% yticks([-1.2:0.3:0.6])
% zticks([0:0.4:0.8])
% xticks([-0.5:0.5:0.5])

view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
set(gca,'FontSize',12)
title('(b)')
ylabel(c,'$\omega_x/|\mathbf{\omega}|$','Interpreter','latex','FontSize',14)
%%
%colormap parula

% f1=sprintf("iso_l2_q2q4_%03d.fig",jcond)
% saveas(h1,f1)
%%

ft=sprintf('../data/velgradfield_dfil_lseQ2_j_%03d.mat',jcond)
md=matfile(ft,'Writable',true)
ftu=sprintf('../data/velgradfield_ufil_lseQ4_j_%03d.mat',jcond);
mu=matfile(ftu,'Writable',true)

h2=figure('OuterPosition',...
    [x1 y1 width height]);
%%
mdlrms=squeeze(rms(md.lambda2tot,[1 2]));
mdlt=-mdlrms(jc);
md.ltot=m2lt;
oz=md.dvdx-md.dudy;%-md.dUdym;
oy=md.dudz-md.dwdx;
ox=md.dwdy-md.dvdz;
omag=sqrt(ox.^2+oy.^2+oz.^2);
mdcos=oz./omag;
subplot(1,2,1)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute(q,[2 1 3]), lt, permute((1e+3)*m.fxQ2,[2 1 3]))

isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(md.lambda2,[2 1 3]), m2lt, permute(mdcos,[2 1 3]))
colormap jet
clim([-1 1])
%colorbar
lightangle(-60,-90)
axis equal

% ylim([-1.501 3.001])
% xlim([-0.601 0.601])
% zlim([0 0.401])
% xticks([-0.6:0.6:0.6])
% yticks([-1.5:0.5:3])
% zticks([0:0.4:0.4])

% ylim([-1.001 2.001])
% xlim([-0.601 0.601])
% zlim([0 0.401])
% xticks([-0.6:0.6:0.6])
% yticks([-1:0.5:2])
% zticks([0:0.4:0.4])

% xlim([-0.801 0.801])
% ylim([-2.5 1.5])
% zlim([0 1.001])
% yticks([-2.5:0.5:1.5])
% zticks([0:0.5:1])
% xticks([-0.8:0.8:0.8])

view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
set(gca,'FontSize',12)
title('(a)')
%%
mulrms=squeeze(rms(mu.lambda2tot,[1 2]));
mult=-mulrms(jc);
mu.ltot=m4lt;
oz=mu.dvdx-mu.dudy;%-mu.dUdym;
oy=mu.dudz-mu.dwdx;
ox=mu.dwdy-mu.dvdz;
omag=sqrt(ox.^2+oy.^2+oz.^2);
mucos=oz./omag;
subplot(1,2,2)
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(mu.lambda2,[2 1 3]), m4lt, permute(mucos,[2 1 3])) 
colormap redblue
clim([-1 1])
cf=colorbar;
lightangle(-60,-90)
axis equal

% ylim([-1.501 3.001])
% xlim([-0.601 0.601])
% zlim([0 0.401])
% xticks([-0.6:0.6:0.6])
% yticks([-1.5:0.5:3])
% zticks([0:0.4:0.4])

% ylim([-1.001 2.001])
% xlim([-0.601 0.601])
% zlim([0 0.401])
% xticks([-0.6:0.6:0.6])
% yticks([-1:0.5:2])
% zticks([0:0.4:0.4])

% xlim([-0.801 0.801])
% ylim([-2.5 1.5])
% zlim([0 1.001])
% yticks([-2.5:0.5:1.5])
% zticks([0:0.5:1])
% xticks([-0.8:0.8:0.8])

view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
set(gca,'FontSize',12)
title('(b)')
colormap redblue
ylabel(cf,'$\omega_x/|\mathbf{\omega}|$','Interpreter','latex','FontSize',14)

% fn=sprintf("iso_l2tot_fil_q2dq4u_j_%03d.fig",jcond)
% saveas(h2,fn)


%%
ft=sprintf('../data/velgradfield_dfil_lseQ2_j_%03d.mat',jcond)
md=matfile(ft,'Writable',true)
ftu=sprintf('../data/velgradfield_ufil_lseQ4_j_%03d.mat',jcond);
mu=matfile(ftu,'Writable',true)

h2=figure('OuterPosition',...
    [x1 y1 width height]);
mdlrms=squeeze(rms(md.lambda2tot,[1 2]));
mdlt=-mdlrms(jc);
md.ltot=m2lt;
oz=md.dvdx-md.dudy-md.dUdym;
oy=md.dudz-md.dwdx;
ox=md.dwdy-md.dvdz;
omag=sqrt(ox.^2+oy.^2+oz.^2);
mdcos=ox./omag;

subplot(1,2,1)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute(q,[2 1 3]), lt, permute((1e+3)*m.fxQ2,[2 1 3]))

isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(md.lambda2tot./ldrms(1,1,111:end),[2 1 3]), -0.05, permute(mdcos,[2 1 3]))
colormap jet
clim([-1 1])
%colorbar
lightangle(-60,-90)
axis equal

ylim([-0.401 0.401])
xlim([-0.401 0.401])
zlim([0 0.401])
xticks([-0.4:0.4:0.4])
yticks([-0.4:0.4:0.4])
zticks([0:0.4:0.4])

% ylim([-0.301 0.301])
% xlim([-0.401 0.401])
% zlim([0 0.401])
% xticks([-0.4:0.4:0.4])
% yticks([-0.3:0.3:0.3])
% zticks([0:0.4:0.4])

% xlim([-0.801 0.801])
% ylim([-0.5 0.5])
% zlim([0.2 0.6])
% yticks([-0.5:0.5:0.5])
% zticks([0.2:0.2:0.6])
% xticks([-0.8:0.4:0.8])

view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
set(gca,'FontSize',12)
title('(a)')
%
mulrms=squeeze(rms(mu.lambda2tot,[1 2]));
mult=-mulrms(jc);
mu.ltot=m4lt;
oz=mu.dvdx-mu.dudy-mu.dUdym;
oy=mu.dudz-mu.dwdx;
ox=mu.dwdy-mu.dvdz;
omag=sqrt(ox.^2+oy.^2+oz.^2);
mucos=ox./omag;
subplot(1,2,2)
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(mu.lambda2tot./lurms(1,1,111:end),[2 1 3]), -0.05, permute(mucos,[2 1 3])) 
colormap redblue
clim([-1 1])
cf=colorbar;
lightangle(-60,-90)
axis equal

ylim([-0.401 0.401])
xlim([-0.401 0.401])
zlim([0 0.401])
xticks([-0.4:0.4:0.4])
yticks([-0.4:0.4:0.4])
zticks([0:0.4:0.4])

% ylim([-0.301 0.301])
% xlim([-0.401 0.401])
% zlim([0 0.401])
% xticks([-0.4:0.4:0.4])
% yticks([-0.3:0.3:0.3])
% zticks([0:0.4:0.4])

% xlim([-0.801 0.801])
% ylim([-0.5 0.5])
% zlim([0.2 0.6])
% yticks([-0.5:0.5:0.5])
% zticks([0.2:0.2:0.6])
% xticks([-0.8:0.4:0.8])

view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
set(gca,'FontSize',12)
title('(b)')
colormap redblue
ylabel(cf,'$\omega_x/|\mathbf{\omega}|$','Interpreter','latex','FontSize',14)

fn=sprintf("iso_l2totrms_fil_q2dq4u_j_%03d.fig",jcond)
saveas(h2,fn)
