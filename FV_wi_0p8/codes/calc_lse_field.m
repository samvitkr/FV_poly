clear
Ny=220;
jcond=156;
jc=Ny-jcond+1;

mp=matfile('../data/mean_profiles.mat')
v= sqrt( 0.5*(mp.vv(jcond,1)+mp.vv(jc,1)) );

fn=sprintf('../data/lse_coeff_B_j_%03d.mat',jcond);
ml=matfile(fn);

fn=sprintf('../data/lse_v_field_B_j_%03d.mat',jcond)
m=matfile(fn,'Writable',true);

m.ulse	=	fftshift(fftshift(v*ml.L11,1),2);
m.vlse	=	fftshift(fftshift(v*ml.L21,1),2);
m.wlse	=	fftshift(fftshift(v*ml.L31,1),2);

m.dudxlse=	fftshift(fftshift(v*ml.L41,1),2);
m.dvdxlse=	fftshift(fftshift(v*ml.L51,1),2);
m.dwdxlse=	fftshift(fftshift(v*ml.L61,1),2);

m.dudylse=	fftshift(fftshift(v*ml.L71,1),2);
m.dvdylse=	fftshift(fftshift(v*ml.L81,1),2);
m.dwdylse=	fftshift(fftshift(v*ml.L91,1),2);

m.dudzlse=	fftshift(fftshift(v*ml.L101,1),2);
m.dvdzlse=	fftshift(fftshift(v*ml.L111,1),2);
m.dwdzlse=	fftshift(fftshift(v*ml.L121,1),2);

m.vozlse=	fftshift(fftshift(v*ml.L131,1),2);
m.woylse=	fftshift(fftshift(v*ml.L141,1),2);

m.fxlse=	fftshift(fftshift(v*ml.L151,1),2);
m.fylse=	fftshift(fftshift(v*ml.L161,1),2);
m.fzlse=	fftshift(fftshift(v*ml.L171,1),2);
