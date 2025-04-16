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
jcond=156;
jc=jcond-Ny/2;
[X,Z,Y]=meshgrid(xp,zp,yp);

plotp=sprintf('../iso_lsevp_delay_j_%03d.fig',jcond)
plotn=sprintf('../iso_lsevn_delay_j_%03d.fig',jcond)


fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
fvgpF=sprintf('../data/lsevp_field_tot_F_j_%03d.mat',jcond)
fvgpB=sprintf('../data/lsevp_field_tot_B_j_%03d.mat',jcond)



%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
m1=matfile(fvgp,'Writable',true);
%m2=matfile(fvgn,'Writable',true);
l1=min(m1.lambda2(1:end,1:end,jc),[],'all');
%l2=min(m2.lambda2(1:end,1:end,jc),[],'all');
val=0.05;


x1=150;
y1=150;
x2=3*450;
y2=350;
hp=figure('OuterPosition',...
    [x1 y1 x2 y2]);
subplot(1,3,2)
m=m1;
clear m1

polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l1,[2 1 3]),val,permute(ox,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,-45)
%camlight('left')
clim([-1 1].*1e-2)
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)


subplot(1,3,3)
m=matfile(fvgpF)
polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
% h1=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l1,[2 1 3]),val,permute(ox,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,-45)
%camlight('left')
clim([-1 1].*1e-2)
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)

subplot(1,3,1)
m=matfile(fvgpB)
polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
% h1=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l1,[2 1 3]),val,permute(ox,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,-45)
%camlight('left')
clim([-1 1].*1e-2)
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)

clear m ox polywork nl

%%
fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
fvgnF=sprintf('../data/lsevn_field_tot_F_j_%03d.mat',jcond)
fvgnB=sprintf('../data/lsevn_field_tot_B_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
m1=matfile(fvgn,'Writable',true);
%m2=matfile(fvgn,'Writable',true);
l1=min(m1.lambda2(1:end,1:end,jc),[],'all');
%l2=min(m2.lambda2(1:end,1:end,jc),[],'all');
val=0.05;


x1=150;
y1=150;
x2=3*450;
y2=350;
hn=figure('OuterPosition',...
    [x1 y1 x2 y2]);
subplot(1,3,2)
m=m1;
clear m1

polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l1,[2 1 3]),val,permute(ox,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,-45)
%camlight('left')
clim([-1 1].*1e-2)
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)


subplot(1,3,3)
m=matfile(fvgnF)
polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
% h1=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l1,[2 1 3]),val,permute(ox,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,-45)
%camlight('left')
clim([-1 1].*1e-2)
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)

subplot(1,3,1)
m=matfile(fvgnB)
polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
% h1=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l1,[2 1 3]),val,permute(ox,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,-45)
%camlight('left')
clim([-1 1].*1e-2)
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)
%%
saveas(hp,plotp)
saveas(hn,plotn)
