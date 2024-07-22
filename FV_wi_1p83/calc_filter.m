ut=278/4667;
dnu=1/278;
Nx=512;
Nz=384;
Ny=220;
Lx=  4*pi;
Lz = 2*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
dkx=kx(2)-kx(1);
dkz=kz(2)-kz(1);
[Kx,Kz]=meshgrid(kx(1:Nx/2),kz(1:Nz/2));
m=matfile('transferfields_mean_visc.mat')

convF=m.vozF-m.woyF;
psc1=convF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
phiconv=real(pscm(1:Nz/2,1:Nx/2,:));
phiconv=0.5*(phiconv+flip(phiconv,3));



psc1=m.vozF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
phivoz=real(pscm(1:Nz/2,1:Nx/2,:));
phivoz=0.5*(phivoz+flip(phivoz,3));


psc1=m.woyF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
phiwoy=real(pscm(1:Nz/2,1:Nx/2,:));
phiwoy=0.5*(phiwoy+flip(phiwoy,3));


m2=matfile('spec2d.mat','Writable',true)
m2.Kx=Kx;
m2.Kz=Kz;
m2.phiconv=phiconv(:,:,1:Ny/2);
m2.phivoz=phivoz(:,:,1:Ny/2);
m2.phiwoy=phiwoy(:,:,1:Ny/2);


