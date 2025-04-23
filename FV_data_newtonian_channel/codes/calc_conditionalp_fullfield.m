close all
clear
%load('ygrid.mat')

Nx=512;
Nz=384;
Ny=220;
Lx=4*pi;
Lz=2*pi;
jcond=171;
jc=Ny-jcond+1;
dx=Lx/Nx;
dz=Lz/Nz;
nf=1;

tstart=0000;
tend=100000;
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
clear mp

un=zeros(Nz,Nx,Ny/2);
vn=zeros(Nz,Nx,Ny/2);
wn=zeros(Nz,Nx,Ny/2);

vozn=zeros(Nz,Nx,Ny/2);
woyn=zeros(Nz,Nx,Ny/2);

s=[Nz Nx];
%vjav=zeros(Nz,Nx);

for time=tstart:tstep:tend
    time
    fvel=sprintf("../data/velfields_%07d.mat",time);
    m=matfile(fvel);

    vj=m.vfield(:,:,jcond);
    vjt=m.vfield(:,:,jc);

    ufieldb=		m.ufield(:,:,Ny/2+1:end);
    ufieldt=flip(	m.ufield(:,:,1:Ny/2),3);
    vfieldb=		m.vfield(:,:,Ny/2+1:end);
    vfieldt=flip(	m.vfield(:,:,1:Ny/2),3);
    wfieldb=		m.wfield(:,:,Ny/2+1:end);
    wfieldt=flip(	m.wfield(:,:,1:Ny/2),3);
    clear m


    ft=sprintf("../data/transferfields_%07d.mat",time);
    mt=matfile(ft);
    vozb=           mt.voz(:,:,Ny/2+1:end);
    vozt=flip(      mt.voz(:,:,1:Ny/2),3);
    woyb=           mt.woy(:,:,Ny/2+1:end);
    woyt=flip(      mt.woy(:,:,1:Ny/2),3);
    clear mt


    %% towards the wall
    %bottom half
        disp('bot')

    [M,I] = max(vj(:));

    [kloc, iloc] = ind2sub(s,I);
    vjc=vj;
    while(abs(M)>abs(vthreshold))
    
        counter=counter+1;
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;

    	vjc=circshift(vjc,[kdelta idelta]);
        vjc(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;

        vjc=circshift(vjc,[-kdelta -idelta]);
        [M,I] = max(vjc(:));
        [kloc, iloc] = ind2sub(s,I);

        un=un		+circshift( ufieldb ,[kdelta idelta]);
        vn=vn		+circshift( vfieldb ,[kdelta idelta]);
        wn=wn		+circshift( wfieldb ,[kdelta idelta]);
    	vozn=vozn	+circshift( vozb ,[kdelta idelta]);
        woyn=woyn	+circshift( woyb ,[kdelta idelta]);

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %top half
    disp('top')
    [M,I] = min(vjt(:));
    [kloc, iloc] = ind2sub(s,I);
    vjc=vjt;

    while(abs(M)>abs(vthreshold))
        counter=counter+1;
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;
    	temp=circshift(vjc,[kdelta idelta]);
        temp(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;
        vjc=circshift(temp,[-kdelta -idelta]);
        [M,I] = min(vjc(:));
        [kloc, iloc] = ind2sub(s,I);

        un=un		+circshift( ufieldt,[kdelta idelta]);
        vn=vn		-circshift( vfieldt,[kdelta idelta]);
        wn=wn		+circshift( wfieldt,[kdelta idelta]);
    	vozn=vozn	+circshift( vozt ,[kdelta idelta]);
        woyn=woyn	+circshift( woyt ,[kdelta idelta]);

    end

end
counter
fc=sprintf("../data/conditionalp_jcond_%03d.mat",jcond);
mc=matfile(fc,'Writable',true);
mc.u=un./counter;
mc.v=vn./counter;
mc.w=wn./counter;
mc.voz=vozn./counter;
mc.woy=woyn./counter;
