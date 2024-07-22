load('ygrid.mat')
nx=512;
nz=384;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx;
zp=lz*[0:nz-1]/nz;
yp=yCheb';
re=4667;
load('dpdx.mat')
%ret=154.4;
%ut=ret/re;
[X,Z,Y]=meshgrid(xp,zp,yp);

time=52000;
ut=mean(ut_ts);
ft=sprintf("transferfields_%07d.mat",time)
mt=matfile(ft)
x1=150;
y1=150;
x2=450;
y2=350;
syz=mt.voz-mt.woy+mt.visc+mt.poly;
syz=(syz./-ut^2);
val=20;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
isosurface(X,Z,(Y+1), abs(syz) ,val, syz)
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
clim([-val val])
colorbar 
%colormap jet
%print(h1,'isotryflp','-dpng');
%saveas(h1,'isotry_0088000.fig')
fp=sprintf("iso_syz_%03d_%07d.fig",val,time)
saveas(h1,fp)
