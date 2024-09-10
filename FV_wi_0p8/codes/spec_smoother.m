function [smoothspec] = spec_smoother(spec,winmax)
% winmax=3;
% md=matfile('../Wi13p5/energ2d.mat');
% j=46;
% spec=md.epF(:,:,j);
[N1,N2]=size(spec);
Nx=2*N1;
Nz=2*N2;
wconv=zeros(winmax+N1,winmax+N2);
pscfull=zeros(Nx,Nz);
pscfull(Nx/2+1:end,Nz/2+1:end)=spec;
pscfull=pscfull+flip(pscfull,1)+flip(pscfull,2)+ flip(flip(pscfull,1),2);
window=winmax;
phishift=zeros(Nx+2*window,Nz+2*window);
phishift(window+1:window+Nx,window+1:window+Nz)=pscfull;

sphi=phishift;
count=1;
for i =1:window

    sphi=sphi+circshift(phishift,i);
    sphi=sphi+circshift(phishift,-i);
    sphi=sphi+circshift(phishift,i,2);
    sphi=sphi+circshift(phishift,-i,2);

    count=count+4;
    for j =1:window
        sphi=sphi+circshift(circshift(phishift,i),j,2);
        sphi=sphi+circshift(circshift(phishift,i),-j,2);
        sphi=sphi+circshift(circshift(phishift,-i),j,2);
        sphi=sphi+circshift(circshift(phishift,-i),j,2);

        count=count+4;

    end
end

sphi=sphi./count;
sphi=0.25*(sphi+flip(sphi,1)+flip(sphi,2)+ flip(flip(sphi,1),2));
        sphiq=sphi(window+Nx/2+1:end,window+Nz/2+1:end);
        [nx, nz]=size(sphiq);
        wconv(1:nx,1:nz)=sphiq;
smoothspec=wconv;
end