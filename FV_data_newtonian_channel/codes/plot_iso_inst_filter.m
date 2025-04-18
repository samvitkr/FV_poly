close all
clear
nx=512;
nz=384;
ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx;%-lx/2;
zp=lz*[0:nz-1]/nz;%-lz/2;
%%
iend=192;
jend=256;
time=40000;
%lt=0.05;
%lt=0.22;vim ca
%lt=-0.001;
lt=-1;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')


mm=matfile('../data/mean_profiles.mat')
dUdy=reshape(mm.dUdy,[1,1,ny]);
dUdy=dUdy(1,1,111:end);
yp=yCheb(111:end)'+1;
[X,Z,Y]=meshgrid(xp,zp,yp);

ft =sprintf('../data/velgrad_transfer_dfil_%07d.mat',time);
%ft =sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond);
%ft=sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond')
m=matfile(ft,'Writable',true)

ftu =sprintf('../data/velgrad_transfer_ufil_%07d.mat',time);
%ftu=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond);
%ftu =sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond);
%ftu=sprintf('../data/velgrad_voz_field_lseQ4_j_156.mat')
mu=matfile(ftu,'Writable',true)

qrmsd=rms(m.lambda2,[1 2]);
qrmsu=rms(mu.lambda2,[1 2]);
l2d=m.lambda2(1:end,1:end,111:end)./qrmsd(1,1,111:end);
l2u=mu.lambda2(1:end,1:end,111:end)./qrmsu(1,1,111:end);

% nld=(m.v).*(m.dvdx-m.dudy)-(m.w).*(m.dudz-m.dwdx);
% nlu=(mu.v).*(mu.dvdx-mu.dudy)-(mu.w).*(mu.dudz-mu.dwdx);
% syzd=nld+m.fx;
% syzu=nlu+mu.fx;
%nld=nld(:,:,111:end);
%nlu=nlu(:,:,111:end);
x1=150;
y1=150;
x2=900;
y2=400;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
%h1=figure;
%subplot(2,1,1)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute(q,[2 1 3]), lt, permute((1e+3)*m.fxQ2,[2 1 3]))
subplot(1,2,1)
isosurface( permute(Z(1:iend,1:jend,1:end),[2 1 3]), permute(X(1:iend,1:jend,1:end),[2 1 3]), permute(Y(1:iend,1:jend,1:end),[2 1 3]),...
permute(l2d(1:iend,1:jend,1:end),[2 1 3]), lt)%, permute(m.v,[2 1 3]))

%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fxQ2,[2 1 3]), -fth)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fxQ2,[2 1 3]), fth)
%clim([-1 1])
colormap jet
%colorbar
lightangle(-45,-90)
axis equal
ylim([0 2*pi])
xlim([0 pi])
zlim([0 1])
zticks([0:0.5:1])
xticks([0:pi:2*pi])
yticks([0:pi:2*pi])
xticklabels({'0','\pi','2\pi'})
yticklabels({'0','\pi','2\pi'})
view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
title('(a)')
%subplot(2,1,2)
%isosurface( permute(mu.Z,[2 1 3]), permute(mu.X,[2 1 3]), permute(mu.Y,[2 1 3]),...
%permute(qu,[2 1 3]), lt, permute((1e+3)*mu.fxQ2,[2 1 3])) 
%%


%figure
subplot(1,2,2)
isosurface( permute(Z(1:iend,1:jend,1:end),[2 1 3]), permute(X(1:iend,1:jend,1:end),[2 1 3]), permute(Y(1:iend,1:jend,1:end),[2 1 3]),...
permute(l2u(1:iend,1:jend,1:end),[2 1 3]), lt)%, permute(mu.v,[2 1 3])) 

%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fyQ2,[2 1 3]), -fth)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fyQ2,[2 1 3]), fth)
%clim([-1 1])
colormap jet
%colorbar
lightangle(-45,-90)
axis equal
ylim([0 2*pi])
xlim([0 pi])
zlim([0 1])
zticks([0:0.5:1])
xticks([0:pi:2*pi])
yticks([0:pi:2*pi])
xticklabels({'0','\pi','2\pi'})
yticklabels({'0','\pi','2\pi'})
view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on
title('(b)')
colormap parula