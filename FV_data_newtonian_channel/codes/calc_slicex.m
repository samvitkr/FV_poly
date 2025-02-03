close all
clear
nx=512;
nz=384;
ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
%jcset=[116 135 187 198 205];
jcset=171;
for jc=1:1
	jcond=jcset(jc)
	%islice=256;
	%lt=0.05;
	%lt=0.22;vim ca
	%lt=-0.005;
	%fth=0.1;
	%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
	%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
	load('../data/ygrid.mat')
	mm=matfile('../data/mean_profiles.mat')
	ut=mm.ret/mm.re;
	dUdy=reshape(mm.dUdy,[1,1,ny]);
	dUdy=dUdy(1,1,111:end);
	
	yp=yCheb(111:end)'+1;
	[Z,Y]=meshgrid(zp,yp);
	
	%ft =sprintf('../data/velgradfield_dfil_lseQ4_j_%03d.mat',jcond);
	ft =sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond);
	%ft=sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond')
	m=matfile(ft,'Writable',true)
	ftu =sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond);
	%ftu=sprintf('../data/velgrad_voz_field_lseQ4_j_156.mat')
	%ftu=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond')
	mu=matfile(ftu,'Writable',true)
	
	for islice=210:280
		fx=sprintf('../data/xslices/velgrad_xslice_lsevpn_i_%03d_j_%03d.mat',islice,jcond)
		mx=matfile(fx,'Writable',true)
		
		ud=  squeeze(m.u(:,islice,:))';	
		vd=  squeeze(m.v(:,islice,:))';
		wd=  squeeze(m.w(:,islice,:))';
		oxd= squeeze(m.dwdy(:,islice,:)-m.dvdz(:,islice,:))';
		ozd= squeeze(m.dvdx(:,islice,:)-m.dudy(:,islice,:)-dUdy)';
		oyd= squeeze(m.dudz(:,islice,:)-m.dwdx(:,islice,:))';
		vozd=squeeze(m.v(:,islice,:).*(m.dvdx(:,islice,:)-m.dudy(:,islice,:)-dUdy))';
		woyd=squeeze(m.w(:,islice,:).*(m.dudz(:,islice,:)-m.dwdx(:,islice,:)))';
		ld=  squeeze(m.lambda2tot(:,islice,:))';
		vozdc=squeeze(m.voz(:,islice,:))';
		woydc=squeeze(m.woy(:,islice,:))';
	
		 uu= squeeze(mu.u(:,islice,:))';	
		vu=  squeeze(mu.v(:,islice,:))';
		wu=  squeeze(mu.w(:,islice,:))';
		oxu= squeeze(mu.dwdy(:,islice,:)-mu.dvdz(:,islice,:))';
		ozu= squeeze(mu.dvdx(:,islice,:)-mu.dudy(:,islice,:)-dUdy)';
		oyu= squeeze(mu.dudz(:,islice,:)-mu.dwdx(:,islice,:))';
		vozu=squeeze(mu.v(:,islice,:).*(mu.dvdx(:,islice,:)-mu.dudy(:,islice,:)-dUdy))';
		woyu=squeeze(mu.w(:,islice,:).*(mu.dudz(:,islice,:)-mu.dwdx(:,islice,:)))';
		lu=  squeeze(mu.lambda2tot(:,islice,:))';
		vozuc=squeeze(mu.voz(:,islice,:))';
	        woyuc=squeeze(mu.woy(:,islice,:))';
	
	
		mx.ud=ud;
		mx.uu=uu;	
		mx.vd=vd;
		mx.vu=vu;
		mx.oxd=oxd;
		mx.oxu=oxu;
		mx.ozd=ozd;
		mx.ozu=ozu;
		mx.wd=wd;
		mx.wu=wu;
		mx.oyd=oyd;
		mx.oyu=oyu;
		mx.vozd=vozd;
		mx.vozu=vozu;
		mx.woyd=woyd;
		mx.woyu=woyu;
		mx.vozdc=vozdc;
	        mx.vozuc=vozuc;
	        mx.woydc=woydc;
	        mx.woyuc=woyuc;
		mx.ld=ld;
		mx.lu=lu;
		mx.Z=Z;
		mx.Y=Y;
	end
end
%%
%  pcolor(Z,Y,woyd)
%  shading flat
%  axis equal
% clim([-1e-5 1e-5])
