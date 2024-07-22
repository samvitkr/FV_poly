load('ygrid.mat')
nx=512;
nz=384;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx;
zp=lz*[0:nz-1]/nz;
yp=yCheb';
re=4667;
ret=150;
ut=ret/re;
[X,Z,Y]=meshgrid(xp,zp,yp);

%mt=matfile('velgrad_transfer_flp_0079000.mat')
%m=matfile('lambdaflp_0079000.mat');
%mt=matfile('polyfields_0050000.mat')

%mt=matfile('transferfields_0079000.mat')
%m=matfile('lambda_0079000.mat');

x1=50;
y1=50;
x2=500;
y2=350;

ts=140000;
istart=140000;
iend=140000;
%for i =istart:iend
	time=ts;
	%time=90000;
	
	ft=sprintf("transferfields_%07d.mat",time)
	mt=matfile(ft);
	f=sprintf("lambda_%07d.mat",time)
	m=matfile(f);
	fn=sprintf("iso_frame_%07d",time)
	
	ml=(mean(mean(m.lambda2,1),2));
	l=m.lambda2;
	lrms=rms(rms(m.lambda2-ml,1),2);
    syz=(mt.voz-mt.woy+mt.visc+mt.poly)./(-ut^2);
	%l=l./lrms;
	h1=figure('OuterPosition',...
	    [x1 y1 x2 y2]);
	%isosurface(X,Z,(Y+1),l,0.05,(mt.voz-mt.woy+mt.visc+mt.poly)./(-ut^2))

    isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y+1,[2 1 3]),...
	permute(l,[2 1 3]), -0.05, permute(syz, [2 1 3]))
	%(mt.voz-mt.woy)./(-ut^2) )
	%pbaspect([4*pi 2*pi 2])
	axis equal
%	clim([-10 10])
%	colorbar 
%	colormap redblue
	
	
	zlim([0.01 1])
	xlabel('$x/H$','interpreter','latex')
	ylabel('$z/H$','interpreter','latex')
	zlabel('$y/H$','interpreter','latex')
	view(-45,45)
	camlight
	lighting gouraud
	lightangle(-45,-90)
%	print(h1,fn,'-dpng');
%	close all
%end
%print(h1,'iso_frame_0079000','-dpng');
saveas(h1,'iso_lambda_n0p05.fig')
