close all
clear
%load('ygrid.mat')

Nx=512;
Nz=384;
Ny=220;
Lx=4*pi;
Lz=2*pi;
jcond=171;
jc=jcond-Ny/2;
dx=Lx/Nx;
dz=Lz/Nz;
nf=1;

tstart=0000;
tend=10000;
%tend=100000;
%tend=0;
tstep=1000;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
itarget=257;
ktarget=193;
wini=round(1/dx);
wink=round(1/dz);
counter=0;
mp=matfile('../data/mean_profiles.mat')
vrms=sqrt(0.5*(mp.vv(jcond,1)+mp.vv(jc,1)));
vthreshold=2*vrms;
un=zeros(Nz,Nx,Ny/2);
vn=zeros(Nz,Nx,Ny/2);
wn=zeros(Nz,Nx,Ny/2);
%vjav=zeros(Nz,Nx);

for time=tstart:tstep:tend
    time
    fvel=sprintf("../data/velfields_%07d.mat",time);
    m=matfile(fvel);
    %     ft=sprintf("../data/transferfields_%07d.mat",time);
    %     mt=matfile(ft);
    %     fg=sprintf("../data/velgrad_%07d.mat",time);
    %     mg=matfile(fg);

    vj=m.vfield(:,:,jcond);
    vjt=m.vfield(:,:,jc);

    ufieldb=m.ufield(:,:,Ny/2+1:end);
    ufieldt=flip(m.ufield(:,:,1:Ny/2));
    vfieldb=m.vfield(:,:,Ny/2+1:end);
    vfieldt=flip(m.vfield(:,:,1:Ny/2));
    wfieldb=m.wfield(:,:,Ny/2+1:end);
    wfieldt=flip(m.wfield(:,:,1:Ny/2));

    clear m
    %uj=m.ufield(:,:,1+Ny/2:end);
    %ujt=flip(m.ufield(:,:,1:Ny/2),3);

    s=size(vj);
    %% towards the wall
%bottom half
    [M,I] = min(vj(:));
    [kloc, iloc] = ind2sub(s,I);
    vjc=vj;
    while(abs(M)>abs(vthreshold))
        counter=counter+1
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;

        un=un+circshift( ufieldb ,[kdelta idelta]);
        vn=vn+circshift( vfieldb ,[kdelta idelta]);
        wn=wn+circshift( wfieldb ,[kdelta idelta]);
        temp=circshift(vjc,[kdelta idelta]);
        temp(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;
        vjc=circshift(temp,[-kdelta -idelta]);
        [M,I] = min(vjc(:));
        [kloc, iloc] = ind2sub(s,I);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%top half   
[M,I] = max(vjt(:));
    [kloc, iloc] = ind2sub(s,I);
    vjc=vjt;
    while(abs(M)>abs(vthreshold))
        counter=counter+1
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;

        un=un+circshift( ufieldt,[kdelta idelta]);
        vn=vn-circshift( vfieldt,[kdelta idelta]);
        wn=wn+circshift( wfieldt,[kdelta idelta]);


        temp=circshift(vjc,[kdelta idelta]);
        temp(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;
        vjc=circshift(temp,[-kdelta -idelta]);
        [M,I] = min(vjc(:));
        [kloc, iloc] = ind2sub(s,I);
    end

end
fc=sprintf("../data/conditionaln_%07d.mat",time);
m=matfile(fc,'Writable',true);
m.u=m.un./counter;
m.v=m.vn./counter;
m.w=m.wn./counter;

