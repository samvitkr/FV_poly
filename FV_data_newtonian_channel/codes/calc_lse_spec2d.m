close all
clear
Nx=512;
Nz=384;
Ny=220;
jcond=156;
jc=jcond-Ny/2;
Lx=  4*pi;
Lz = 2*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
dkx=kx(2)-kx(1);
dkz=kz(2)-kz(1);
[Kx,Kz]=meshgrid(kx(1:Nx/2),kz(1:Nz/2));
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
%ft =sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond);
ft =sprintf('../data/velgradfx_voz_field_lseQ2_j_%03d.mat',jcond);
m2=matfile(ft)
%ftu=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond);
ftu =sprintf('../data/velgradfx_voz_field_lseQ4_j_%03d.mat',jcond);
m4=matfile(ftu)


nl2F=conj(fft2(m2.v)).*fft2(m2.dvdx-m2.dudz)-conj(fft2(m2.w)).*fft2(m2.dudz-m2.dwdx);
nl4F=conj(fft2(m4.v)).*fft2(m4.dvdx-m4.dudz)-conj(fft2(m4.w)).*fft2(m4.dudz-m4.dwdx);
nl2F=nl2F./(Nx^2*Nz^2);
nl4F=nl4F./(Nx^2*Nz^2);

psc1=nl2F;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
phinl2=real(pscm(1:Nz/2,1:Nx/2,:));

psc1=nl4F;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
phinl4=real(pscm(1:Nz/2,1:Nx/2,:));


fpu=sprintf('../data/velgrad_lse_voz_spec_j_%03d.mat',jcond);
mpu=matfile(fpu,'Writable',true)
mpu.phinl2=phinl2;
mpu.phinl4=phinl4;
mpu.Kx=Kx;
mpu.Kz=Kz;
