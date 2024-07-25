load('ygrid.mat')
nx=512;
nz=384;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx;
zp=lz*[0:nz-1]/nz;
yp=yCheb';
re=4667;
ret=232;
ut=ret/re;
[X,Z,Y]=meshgrid(xp,zp,yp);
mt=matfile('velgrad_transfer_flp_0070000.mat')
m=matfile('lambdaflp_0070000.mat');
%mt=matfile('transferfields_0040000.mat')
%m=matfile('lambda_0040000.mat');
ml=(mean(mean(m.lambda2,1),2));
l=m.lambda2;
lrms=rms(rms(m.lambda2-ml,1),2);
x1=150;
y1=150;
x2=450;
y2=350;
val=-0.02;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
isosurface(X,Z,Y+1,m.lambda2./lrms,val,(mt.voz-mt.woy)./(-ut^2))
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
clim([-10 10])
colorbar 
colormap jet
%print(h1,'isotryflp','-dpng');
saveas(h1,'iso_lambda_flp_70000.fig')
