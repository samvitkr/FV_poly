clear
Ny=220;
jcond=156;
jc=Ny-jcond+1;

mp=matfile('../data/mean_profiles.mat')
U=0.5*(mp.Um(Ny/2+1:end,1)+flip( mp.Um(1:Ny/2,1) ,1));
voz=0.5*(mp.vozm(Ny/2+1:end,1)+flip( mp.vozm(1:Ny/2,1) ,1));
woy=0.5*(mp.woym(Ny/2+1:end,1)+flip( mp.woym(1:Ny/2,1) ,1));
polyfx=0.5*(mp.polym(Ny/2+1:end,1)+flip( mp.polym(1:Ny/2,1) ,1));
dUdy=0.5*(mp.dUdy(Ny/2+1:end,1)-flip( mp.dUdy(1:Ny/2,1) ,1));

U=reshape(U,[1 1 Ny/2]);
voz=reshape(voz,[1 1 Ny/2]);
woy=reshape(woy,[1 1 Ny/2]);
polyfx=reshape(polyfx,[1 1 Ny/2]);
dUdy=reshape(dUdy,[1 1 Ny/2]);

clear mp

%fn=sprintf('../data/lse_v_field_j_%03d.mat',jcond)
%m=matfile(fn,'Writable',true)

%fnp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%mp=matfile(fnp,'Writable',true);

%fn=sprintf('../data/lse_v_field_F_j_%03d.mat',jcond)
%m=matfile(fn,'Writable',true)

%fnp=sprintf('../data/lsevp_field_tot_F_j_%03d.mat',jcond)
%mp=matfile(fnp,'Writable',true);

fn=sprintf('../data/lse_v_field_B_j_%03d.mat',jcond)
m=matfile(fn,'Writable',true)

fnp=sprintf('../data/lsevp_field_tot_B_j_%03d.mat',jcond)
mp=matfile(fnp,'Writable',true);
mp.u=U+	m.ulse;
mp.v=	m.vlse;
mp.w=	m.wlse;

mp.dudx=m.dudxlse;
mp.dvdx=m.dvdxlse;
mp.dwdx=m.dwdxlse;

mp.dudy=dUdy+	m.dudylse;
mp.dvdy=	m.dvdylse;
mp.dwdy=	m.dwdylse;

mp.dudz=m.dudzlse;
mp.dvdz=m.dvdzlse;
mp.dwdz=m.dwdzlse;

mp.voz=voz+	m.vozlse;
mp.woy=woy+	m.woylse;

mp.fx=polyfx+	m.fxlse;
mp.fy=		m.fylse;
mp.fz=		m.fzlse;

clear mp

%fnn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
%mn=matfile(fnn,'Writable',true);


%fnn=sprintf('../data/lsevn_field_tot_F_j_%03d.mat',jcond)
%mn=matfile(fnn,'Writable',true);

fnn=sprintf('../data/lsevn_field_tot_B_j_%03d.mat',jcond)
mn=matfile(fnn,'Writable',true);

mn.u=U- m.ulse;
mn.v= - m.vlse;
mn.w= - m.wlse;

mn.dudx=-m.dudxlse;
mn.dvdx=-m.dvdxlse;
mn.dwdx=-m.dwdxlse;

mn.dudy=dUdy-   m.dudylse;
mn.dvdy=    -   m.dvdylse;
mn.dwdy=    -   m.dwdylse;

mn.dudz=-m.dudzlse;
mn.dvdz=-m.dvdzlse;
mn.dwdz=-m.dwdzlse;

mn.voz=voz-     m.vozlse;
mn.woy=woy-     m.woylse;

mn.fx=polyfx-     m.fxlse;
mn.fy=    -     m.fylse;
mn.fz=    -     m.fzlse;

clear mn

clear m
