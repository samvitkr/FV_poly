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
jcset=180;
for jc=1:1
	jcond=jcset(jc)
	load('../data/ygrid.mat')
	mm=matfile('../data/mean_profiles.mat')
	ut=mm.ret/mm.re;
	yp=yCheb(111:end)'+1;
	[Z,Y]=meshgrid(zp,yp);
	
	ft=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
	m=matfile(ft,'Writable',true)
	ftu=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
	mu=matfile(ftu,'Writable',true)
	
	islice=257;
		fx=sprintf('../data/xslices/vel_lsevpn_i_%03d_j_%03d.mat',islice,jcond)
		mx=matfile(fx,'Writable',true)

		ud=  squeeze(m.u(:,islice,:))';	
		vd=  squeeze(m.v(:,islice,:))';
		wd=  squeeze(m.w(:,islice,:))';
		oxd= squeeze(m.dwdy(:,islice,:)-m.dvdz(:,islice,:))';
		ozd= squeeze(m.dvdx(:,islice,:)-m.dudy(:,islice,:))';
		oyd= squeeze(m.dudz(:,islice,:)-m.dwdx(:,islice,:))';
		vozd=squeeze(m.v(:,islice,:).*(m.dvdx(:,islice,:)-m.dudy(:,islice,:)))';
		woyd=squeeze(m.w(:,islice,:).*(m.dudz(:,islice,:)-m.dwdx(:,islice,:)))';
	
		vozdc=squeeze(m.voz(:,islice,:))';
		woydc=squeeze(m.woy(:,islice,:))';
		fxd=squeeze(m.fx(:,islice,:))';
		fyd=squeeze(m.fy(:,islice,:))';
		fzd=squeeze(m.fz(:,islice,:))';

		uu= squeeze(mu.u(:,islice,:))';	
		vu= squeeze(mu.v(:,islice,:))';
		wu= squeeze(mu.w(:,islice,:))';
		oxu=squeeze(mu.dwdy(:,islice,:)-mu.dvdz(:,islice,:))';
		ozu=squeeze(mu.dvdx(:,islice,:)-mu.dudy(:,islice,:))';
		oyu=squeeze(mu.dudz(:,islice,:)-mu.dwdx(:,islice,:))';
		vozu=squeeze(mu.v(:,islice,:).*(mu.dvdx(:,islice,:)-mu.dudy(:,islice,:)))';
		woyu=squeeze(mu.w(:,islice,:).*(mu.dudz(:,islice,:)-mu.dwdx(:,islice,:)))';
		vozuc=squeeze(mu.voz(:,islice,:))';
	        woyuc=squeeze(mu.woy(:,islice,:))';
		fxu=squeeze(mu.fx(:,islice,:))';
                fyu=squeeze(mu.fy(:,islice,:))';
                fzu=squeeze(mu.fz(:,islice,:))';
	
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

		mx.fxd=fxd;
		mx.fxu=fxu;
		mx.fyd=fyd;
                mx.fyu=fyu;
		mx.fzd=fzd;
                mx.fzu=fzu;

		mx.Z=Z;
		mx.Y=Y;
end
%%
