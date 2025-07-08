close all
clear
load('../data/ygrid.mat')
load('../data/mean_profiles.mat')
ut=ret/re;

nx=512;
nz=384;
Ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
yp=yCheb(Ny/2+1:end)'+1;

jcond=171;
jc=jcond-Ny/2;
%fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
fvgp=sprintf("../data/conditionalp_jcond_%03d.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_%03d.mat",jcond);
m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
[nzz nxx nyy]=size(m1.lambda2);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;
itarget=257;
ktarget=193;
xp=xp(itarget-wxx:itarget+wxx);
zp=zp(ktarget-wzz:ktarget+wzz);
l1=min(m1.lambda2(1:end,1:end,jc),[],'all');
l2=min(m2.lambda2(1:end,1:end,jc),[],'all');
val=-3*ut^2/yp(jc)^2;
[X,Z,Y]=meshgrid(xp,zp,yp);
m1.torx=Y.*m1.polyz-Z.*m1.polyy;
m1.tory=Z.*m1.polyx-X.*m1.polyz;
m1.torz=X.*m1.polyy-Y.*m1.polyx;
m1.ox=m1.dwdy-m1.dvdz;
m1.oy=m1.dudz-m1.dwdx;
m1.oz=m1.dvdx-m1.dudy;
m1.cos_tor_vor = (m1.torx.*m1.ox + m1.tory.*m1.oy +m1.torz.*m1.oz)...
		./sqrt( (m1.torx.^2 + m1.tory.^2 + m1.torz.^2).*...
			(m1.ox.^2 + m1.oy.^2 + m1.oz.^2));

m2.torx=Y.*m2.polyz-Z.*m2.polyy;
m2.tory=Z.*m2.polyx-X.*m2.polyz;
m2.torz=X.*m2.polyy-Y.*m2.polyx;
m2.ox=m2.dwdy-m2.dvdz;
m2.oy=m2.dudz-m2.dwdx;
m2.oz=m2.dvdx-m2.dudy;
m2.cos_tor_vor = (m2.torx.*m2.ox + m2.tory.*m2.oy +m2.torz.*m2.oz)...
                ./sqrt( (m2.torx.^2 + m2.tory.^2 + m2.torz.^2).*...
                        (m2.ox.^2 + m2.oy.^2 + m2.oz.^2));
%
%x1=150;
%y1=-150;
%x2=2*450;
%y2=350;
%h1=figure('OuterPosition',...
%    [x1 y1 x2 y2]);
%subplot(1,2,1)
%m=m1;
%%polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
%%nl=m.woy-m.voz;
%ox=m.dwdy-m.dvdz;
%subplot(1,2,1)
%hold on
%isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2,[2 1 3]),val,permute(ox,[2 1 3]))
%scatter3(0,0,yp(jc),80,'green','filled')
%hold off
%%(mt.voz-mt.woy)./(-ut^2) )
%%pbaspect([4*pi 2*pi 2])
%axis equal
%axis tight
%shading flat
%lightangle(45,-45)
%%camlight('left')
%clim([-1 1])
%colorbar 
%colormap redblue
%%print(h1,'isotryflp','-dpng');
%%saveas(h1,'iso_lambda_flp_70000.fig')
%xlabel('z')
%ylabel('x')
%zlabel('y')
%view(45,45)
%subplot(1,2,2)
%
%m=m2;
%%polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
%%nl=m.woy-m.voz;
%ox=m.dwdy-m.dvdz;
%% h1=figure('OuterPosition',...
%%     [x1 y1 x2 y2]);
%hold on
%isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2,[2 1 3]),val,permute(ox,[2 1 3]))
%scatter3(0,0,yp(jc),80,'green','filled')
%hold off
%%(mt.voz-mt.woy)./(-ut^2) )
%%pbaspect([4*pi 2*pi 2])
%axis equal
%axis tight
%shading flat
%lightangle(45,-45)
%%camlight('left')
%clim([-1 1])
%colorbar 
%colormap redblue
%%print(h1,'isotryflp','-dpng');
%%saveas(h1,'iso_lambda_flp_70000.fig')
%xlabel('z')
%ylabel('x')
%zlabel('y')
%view(45,45)
%
%saveas(h1,'iso_cond_0p8.fig')
