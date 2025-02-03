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

for time=0000:1000:120000
	ut=ut_ts(time+1);
	ft=sprintf("transferfields_%07d.mat",time)
	mt=matfile(ft)
	x1=150;
	y1=150;
	x2=450;
	y2=350;
	syz=mt.voz-mt.woy+mt.visc+mt.poly;

	f=sprintf("lambda_%07d.mat",time)
	m=matfile(f);

	ft1=sprintf("transferfieldsav_yz_%07d.mat",time)
	mt1=matfile(ft1,'Writable',true)
	mt1.voz=squeeze(mean( mt.voz,2  ));
	mt1.woy=squeeze(mean( mt.woy,2  ));
	mt1.visc=squeeze(mean( mt.visc,2  ));
	mt1.poly=squeeze(mean( mt.poly,2  ));
	mt1.syz=squeeze(mean( syz,2  ));
	mt1.l=squeeze(mean( m.lambda2,2  ));
	mt1.Q=squeeze(mean( m.Q,2  ));
end
%syz=(syz./-ut^2);
%val=10;
%h1=figure('OuterPosition',...
%    [x1 y1 x2 y2]);
%isosurface(X,Z,(Y+1), abs(syz) ,val, syz)
% %(mt.voz-mt.woy)./(-ut^2) )
% %pbaspect([4*pi 2*pi 2])
%axis equal
%clim([-val val])
%colorbar 
% %colormap jet
% %print(h1,'isotryflp','-dpng');
% %saveas(h1,'isotry_0088000.fig')
%fp=sprintf("iso_syz_%03d_%07d.fig",val,time)
%saveas(h1,fp)
