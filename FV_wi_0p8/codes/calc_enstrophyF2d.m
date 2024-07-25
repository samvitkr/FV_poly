ret=259;
re=4667;
ut=ret/re;
dnu=1/ret;
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
m=matfile('enstrophy_transferF.mat')


psc1=m.pioxF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
pix=real(pscm(1:Nz/2,1:Nx/2,:));
pix=0.5*(pix+flip(pix,3));


psc1=m.pioyF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
piy=real(pscm(1:Nz/2,1:Nx/2,:));
piy=0.5*(piy+flip(piy,3));

psc1=m.piozF;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
piz=real(pscm(1:Nz/2,1:Nx/2,:));
piz=0.5*(piz+flip(piz,3));



m2=matfile('enstrophy2d.mat','Writable',true)
m2.Kx=Kx;
m2.Kz=Kz;
%m2.epF=phiconv(:,:,1:Ny/2);
m2.pioxF=pix(:,:,1:Ny/2);
m2.pioyF=piy(:,:,1:Ny/2);
m2.piozF=piz(:,:,1:Ny/2);

m2.piox=squeeze(mean(m.piox,[1 2]));
m2.pioy=squeeze(mean(m.pioy,[1 2]));
m2.pioz=squeeze(mean(m.pioz,[1 2]));
%m2.dis=m.dis;
%m2.phivoz=phivoz(:,:,1:Ny/2);
%m2.phiwoy=phiwoy(:,:,1:Ny/2);


