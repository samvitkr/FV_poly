re=4667;
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
m=matfile('energy_transferF.mat')

convF=m.epF;
psc1=convF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
phiconv=real(pscm(1:Nz/2,1:Nx/2,:));
phiconv=0.5*(phiconv+flip(phiconv,3));



psc1=m.pixF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
pix=real(pscm(1:Nz/2,1:Nx/2,:));
pix=0.5*(pix+flip(pix,3));


psc1=m.piyF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
piy=real(pscm(1:Nz/2,1:Nx/2,:));
piy=0.5*(piy+flip(piy,3));

psc1=m.pizF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
piz=real(pscm(1:Nz/2,1:Nx/2,:));
piz=0.5*(piz+flip(piz,3));



m2=matfile('energ2d.mat','Writable',true)
m2.Kx=Kx;
m2.Kz=Kz;
m2.epF=phiconv(:,:,1:Ny/2);
m2.pixF=pix(:,:,1:Ny/2);
m2.piyF=piy(:,:,1:Ny/2);
m2.pizF=piz(:,:,1:Ny/2);
m2.pix=squeeze(mean(m.pix,[1,2]));
m2.piy=squeeze(mean(m.piy,[1,2]));
m2.piz=squeeze(mean(m.piz,[1,2]));

m2.dis=m.dis;
%m2.phivoz=phivoz(:,:,1:Ny/2);
%m2.phiwoy=phiwoy(:,:,1:Ny/2);


