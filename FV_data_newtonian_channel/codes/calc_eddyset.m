close all
clear
nx=512;
nz=384;
ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
%jcset=[116 135 187 198 205]
jcset=171;
%jcond=184;
%jc=jcond-110;
%lt=0.05;
%lt=0.22;vim ca
%lt=-0.005;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
mm=matfile('../data/mean_profiles.mat')
dUdy=reshape(mm.dUdy,[1,1,ny]);
dUdy=dUdy(1,1,111:end);
yp=yCheb(111:end)'+1;
[X,Z,Y]=meshgrid(xp,zp,yp);
lt=-0.005;
%ft =sprintf('../data/velgradfield_dfil_lseQ4_j_%03d.mat',jcond);


for  jj=1:1
	jcond=jcset(jj);
	jc=jcond-110;

	ft =sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond);
	%ft=sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond')
	m=matfile(ft,'Writable',true)
	
	%ftu=sprintf('../data/velgradfield_ufil_lseQ4_j_%03d.mat',jcond);
	%ftu=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond);
	ftu =sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond);
	%ftu=sprintf('../data/velgrad_voz_field_lseQ4_j_156.mat')
	mu=matfile(ftu,'Writable',true)
	
	%qrmsd=rms(m.Q,[1 2]);
	%qrmsu=rms(mu.Q,[1 2]);
	
	nld=(m.v).*(m.dvdx-m.dudy)-(m.w).*(m.dudz-m.dwdx);
	nlu=(mu.v).*(mu.dvdx-mu.dudy)-(mu.w).*(mu.dudz-mu.dwdx);
	syzd=nld+m.fx;
	syzu=nlu+mu.fx;
	%nld=nld(:,:,111:end);
	%nlu=nlu(:,:,111:end);
	islice=[206:306];
	kslice=[142:242];
	ud=  squeeze(m.u(kslice,islice,:));	
	vd=  squeeze(m.v(kslice,islice,:));
	wd=  squeeze(m.w(kslice,islice,:));
	oxd= squeeze(m.dwdy(kslice,islice,:)-m.dvdz(kslice,islice,:));
	ozd= squeeze(m.dvdx(kslice,islice,:)-m.dudy(kslice,islice,:)-dUdy);
	oyd= squeeze(m.dudz(kslice,islice,:)-m.dwdx(kslice,islice,:));
	vozd=squeeze(m.v(kslice,islice,:).*(m.dvdx(kslice,islice,:)-m.dudy(kslice,islice,:)-dUdy));
	woyd=squeeze(m.w(kslice,islice,:).*(m.dudz(kslice,islice,:)-m.dwdx(kslice,islice,:)));
	ld=  squeeze(m.lambda2tot(kslice,islice,:));
	vozdc=squeeze(m.voz(kslice,islice,:));
	woydc=squeeze(m.woy(kslice,islice,:));
	
	
	uu=  squeeze(mu.u(kslice,islice,:));
	vu=  squeeze(mu.v(kslice,islice,:));
	wu=  squeeze(mu.w(kslice,islice,:));
	oxu= squeeze(mu.dwdy(kslice,islice,:)-mu.dvdz(kslice,islice,:));
	ozu= squeeze(mu.dvdx(kslice,islice,:)-mu.dudy(kslice,islice,:)-dUdy);
	oyu= squeeze(mu.dudz(kslice,islice,:)-mu.dwdx(kslice,islice,:));
	vozu=squeeze(mu.v(kslice,islice,:).*(mu.dvdx(kslice,islice,:)-mu.dudy(kslice,islice,:)-dUdy));
	woyu=squeeze(mu.w(kslice,islice,:).*(mu.dudz(kslice,islice,:)-mu.dwdx(kslice,islice,:)));
	lu=  squeeze(mu.lambda2tot(kslice,islice,:));
	vozuc=squeeze(mu.voz(kslice,islice,:));
	woyuc=squeeze(mu.woy(kslice,islice,:));
	
	x=m.X(kslice,islice,:);
	y=m.Y(kslice,islice,:);
	z=m.Z(kslice,islice,:);
	[startX,startY]=meshgrid(0.4,0.05:0.05:0.8);
	
	msn=sprintf('../data/eddyset_d_j_%03d.mat',jcond)
	ms=matfile(msn,'Writable',true)
	ms.ud=permute(ud,[2 1 3]);
	ms.vd=permute(vd,[2 1 3]);
	ms.wd=permute(wd,[2 1 3]);
	ms.oxd=permute(oxd,[2 1 3]);
	ms.oyd=permute(oyd,[2 1 3]);
	ms.ozd=permute(ozd,[2 1 3]);
	ms.ld=permute(ld,[2 1 3]);
	ms.vozdc=permute(vozdc,[2 1 3]);
	ms.woydc=permute(woydc,[2 1 3]);
	
	ms.uu=permute(uu,[2 1 3]);
	ms.vu=permute(vu,[2 1 3]);
	ms.wu=permute(wu,[2 1 3]);
	ms.oxu=permute(oxu,[2 1 3]);
	ms.oyu=permute(oyu,[2 1 3]);
	ms.ozu=permute(ozu,[2 1 3]);
	ms.lu=permute(lu,[2 1 3]);
	ms.vozuc=permute(vozuc,[2 1 3]);
	ms.woyuc=permute(woyuc,[2 1 3]);
	
	ms.x=permute(x,[2 1 3]);
	ms.y=permute(y,[2 1 3]);
	ms.z=permute(z,[2 1 3]);
end

%[sz,sx,sy]=meshgrid(0.4:0.1:0.5,0.4:0.1:0.5,0.2:0.1:0.5);
%vertsv=stream2(squeeze(z(:,50,:)),squeeze(y(:,50,:)),squeeze(ozd(:,50,:)),squeeze(oyd(:,50,:)),startX,startY);
%streamline(vertsv);
%verts=stream3(z,x,y,ozd,oxd,oyd,0.4,0,0.2);
%lineobj=streamline(verts);
%view(3)

% figure
% pcolor(squeeze(ld(:,55,:)))
% shading flat
% clim([-0.005 0.005])
%%
% x1=150;
% y1=150;
% x2=1000;
% y2=350;
% h1=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
% subplot(1,2,1)
% %isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% %permute(q,[2 1 3]), lt, permute((1e+3)*m.fxQ2,[2 1 3]))
% 
% isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
% permute(m.lambda2tot,[2 1 3]), lt, permute(m.v,[2 1 3]))
% 
% %isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% %permute((1e+3)*m.fxQ2,[2 1 3]), -fth)
% %isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% %permute((1e+3)*m.fxQ2,[2 1 3]), fth)
% %clim([-1 1])
% colormap jet
% colorbar
% lightangle(-45,-90)
% axis equal
%  ylim([-1 0.6])
%  xlim([-0.5 0.5])
%  zlim([0 0.8])
% view(45,45)
% xlabel('z')
% ylabel('x')
% zlabel('y')
% grid on
% 
% subplot(1,2,2)
% %isosurface( permute(mu.Z,[2 1 3]), permute(mu.X,[2 1 3]), permute(mu.Y,[2 1 3]),...
% %permute(qu,[2 1 3]), lt, permute((1e+3)*mu.fxQ2,[2 1 3])) 
% 
% isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
% permute(mu.lambda2tot,[2 1 3]), lt, permute(mu.v,[2 1 3])) 
% 
% %isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% %permute((1e+3)*m.fyQ2,[2 1 3]), -fth)
% %isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% %permute((1e+3)*m.fyQ2,[2 1 3]), fth)
% %clim([-1 1])
% colormap jet
% colorbar
% lightangle(-45,-90)
% axis equal
%  ylim([-1 0.6])
%  xlim([-0.5 0.5])
%  zlim([0 0.8])
% view(45,45)
% xlabel('z')
% ylabel('x')
% zlabel('y')
% grid on
% %%

