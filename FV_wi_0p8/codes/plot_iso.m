close all
clear
load('../data/ygrid.mat')
nx=512;
nz=384;
Ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
yp=yCheb(Ny/2+1:end)'+1;
jcond=180;
jc=jcond-Ny/2;
fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)

m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
l1=min(m1.lambda2(1:end,1:end,jc),[],'all');
l2=min(m2.lambda2(1:end,1:end,jc),[],'all');
val=0.05;
[X,Z,Y]=meshgrid(xp,zp,yp);
% mt=matfile('velgrad_transfer_flp_0070000.mat')
% m=matfile('lambdaflp_0070000.mat');
%mt=matfile('transferfields_0040000.mat')
%m=matfile('lambda_0040000.mat');
% ml=(mean(mean(m.lambda2,1),2));
% l=m.lambda2;
% lrms=rms(rms(m.lambda2-ml,1),2);

x1=150;
y1=150;
x2=2*450;
y2=350;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
subplot(1,2,1)
m=m1;
polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
subplot(1,2,1)
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l2,[2 1 3]),val,permute(polywork,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,-45)
%camlight('left')
clim([-1 1].*1e-4)
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)
subplot(1,2,2)

m=m2;
polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
% h1=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l2,[2 1 3]),val,permute(polywork,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,-45)
%camlight('left')
clim([-1 1].*1e-4)
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)
