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

tstart=40000;
tend=130000;
%tend=100000;
%tend=0;
tstep=1000;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
itarget=Nx/2+1;
ktarget=Nz/2+1;
xbox=1;
zbox=0.8;
wini=round(xbox/dx);
wink=round(zbox/dz);
winiav=round(0.5*wini);
winkav=round(0.5*wink);
nzav=2*winkav+1;
nxav=2*winiav+1;
event_location=[];

counter=0;
mp=matfile('../data/mean_profiles.mat')
vrms=sqrt(0.5*(mp.vv(jcond,1)+mp.vv(jc,1)));
vthreshold=vrms;
clear mp

un=zeros(nzav,nxav,Ny/2);
vn=zeros(nzav,nxav,Ny/2);
wn=zeros(nzav,nxav,Ny/2);

dudxn=zeros(nzav,nxav,Ny/2);
dvdxn=zeros(nzav,nxav,Ny/2);
dwdxn=zeros(nzav,nxav,Ny/2);

dudyn=zeros(nzav,nxav,Ny/2);
dvdyn=zeros(nzav,nxav,Ny/2);
dwdyn=zeros(nzav,nxav,Ny/2);

dudzn=zeros(nzav,nxav,Ny/2);
dvdzn=zeros(nzav,nxav,Ny/2);
dwdzn=zeros(nzav,nxav,Ny/2);

vozn=zeros(nzav,nxav,Ny/2);
woyn=zeros(nzav,nxav,Ny/2);

polyxn=zeros(nzav,nxav,Ny/2);
polyyn=zeros(nzav,nxav,Ny/2);
polyzn=zeros(nzav,nxav,Ny/2);

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


    fvelg=sprintf("../data/velgrad_%07d.mat",time);
    mg=matfile(fvelg)
    dudxb=            mg.dudx(:,:,Ny/2+1:end);
    dudxt=flip(       mg.dudx(:,:,1:Ny/2),3);
    dvdxb=            mg.dvdx(:,:,Ny/2+1:end);
    dvdxt=flip(       mg.dvdx(:,:,1:Ny/2),3);
    dwdxb=            mg.dwdx(:,:,Ny/2+1:end);
    dwdxt=flip(       mg.dwdx(:,:,1:Ny/2),3);

    dudyb=            mg.dudy(:,:,Ny/2+1:end);
    dudyt=flip(       mg.dudy(:,:,1:Ny/2),3);
    dvdyb=            mg.dvdy(:,:,Ny/2+1:end);
    dvdyt=flip(       mg.dvdy(:,:,1:Ny/2),3);
    dwdyb=            mg.dwdy(:,:,Ny/2+1:end);
    dwdyt=flip(       mg.dwdy(:,:,1:Ny/2),3);

    dudzb=            mg.dudz(:,:,Ny/2+1:end);
    dudzt=flip(       mg.dudz(:,:,1:Ny/2),3);
    dvdzb=            mg.dvdz(:,:,Ny/2+1:end);
    dvdzt=flip(       mg.dvdz(:,:,1:Ny/2),3);
    dwdzb=            mg.dwdz(:,:,Ny/2+1:end);
    dwdzt=flip(       mg.dwdz(:,:,1:Ny/2),3);
    clear mg


    ft=sprintf("../data/transferfields_%07d.mat",time);
    mt=matfile(ft);
    vozb=           mt.voz(:,:,Ny/2+1:end);
    vozt=flip(      mt.voz(:,:,1:Ny/2),3);
    woyb=           mt.woy(:,:,Ny/2+1:end);
    woyt=flip(      mt.woy(:,:,1:Ny/2),3);

    polyxb=         mt.polyx(:,:,Ny/2+1:end);
    polyxt=flip(    mt.polyx(:,:,1:Ny/2),3);
    polyyb=         mt.polyy(:,:,Ny/2+1:end);
    polyyt=flip(    mt.polyy(:,:,1:Ny/2),3);
    polyzb=         mt.polyz(:,:,Ny/2+1:end);
    polyzt=flip(    mt.polyz(:,:,1:Ny/2),3);

    clear mt



    %% towards the wall
    %bottom half
    disp('bot')

    [M,I] = min(vj(:));

    [kloc, iloc] = ind2sub(s,I);
    vjc=vj;
    while(abs(M)>abs(vthreshold))
	event_location=[event_location;kloc iloc jcond time];

        counter=counter+1;
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;

        vjc=circshift(vjc,[kdelta idelta]);
        vjc(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;

        vjc=circshift(vjc,[-kdelta -idelta]);
        [M,I] = min(vjc(:));
        [kloc, iloc] = ind2sub(s,I);

        ufieldb=circshift( ufieldb ,[kdelta idelta]);
        vfieldb=circshift( vfieldb ,[kdelta idelta]);
        wfieldb=circshift( wfieldb ,[kdelta idelta]);

        dudxb=circshift( dudxb ,[kdelta idelta]);
        dvdxb=circshift( dvdxb ,[kdelta idelta]);
        dwdxb=circshift( dwdxb ,[kdelta idelta]);

        dudyb=circshift( dudyb ,[kdelta idelta]);
        dvdyb=circshift( dvdyb ,[kdelta idelta]);
        dwdyb=circshift( dwdyb ,[kdelta idelta]);

        dudzb=circshift( dudzb ,[kdelta idelta]);
        dvdzb=circshift( dvdzb ,[kdelta idelta]);
        dwdzb=circshift( dwdzb ,[kdelta idelta]);

        vozb=circshift( vozb ,[kdelta idelta]);
        woyb=circshift( woyb ,[kdelta idelta]);

	polyxb=circshift( polyxb ,[kdelta idelta]);
	polyyb=circshift( polyyb ,[kdelta idelta]);
	polyzb=circshift( polyzb ,[kdelta idelta]);

        un=	un      +ufieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        vn=	vn      +vfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        wn=	wn      +wfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudxn=dudxn     +dudxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdxn=dvdxn     +dvdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdxn=dwdxn     +dwdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudyn=dudyn     +dudyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdyn=dvdyn     +dvdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdyn=dwdyn     +dwdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudzn=dudzn     +dudzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdzn=dvdzn     +dvdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdzn=dwdzn     +dwdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        vozn=	vozn	+vozb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        woyn=	woyn	+woyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

	polyxn=  polyxn   +polyxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
	polyyn=  polyyn   +polyyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
	polyzn=  polyzn   +polyzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        ufieldb=circshift( ufieldb ,-[kdelta idelta]);
        vfieldb=circshift( vfieldb ,-[kdelta idelta]);
        wfieldb=circshift( wfieldb ,-[kdelta idelta]);
        dudxb=circshift( dudxb ,-[kdelta idelta]);
        dvdxb=circshift( dvdxb ,-[kdelta idelta]);
        dwdxb=circshift( dwdxb ,-[kdelta idelta]);
        dudyb=circshift( dudyb ,-[kdelta idelta]);
        dvdyb=circshift( dvdyb ,-[kdelta idelta]);
        dwdyb=circshift( dwdyb ,-[kdelta idelta]);
        dudzb=circshift( dudzb ,-[kdelta idelta]);
        dvdzb=circshift( dvdzb ,-[kdelta idelta]);
        dwdzb=circshift( dwdzb ,-[kdelta idelta]);
        vozb=circshift( vozb ,-[kdelta idelta]);
        woyb=circshift( woyb ,-[kdelta idelta]);

	polyxb=circshift( polyxb ,-[kdelta idelta]);
	polyyb=circshift( polyyb ,-[kdelta idelta]);
	polyzb=circshift( polyzb ,-[kdelta idelta]);

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    top half
    disp('top')
    [M,I] = max(vjt(:));
    [kloc, iloc] = ind2sub(s,I);
    vjc=vjt;

    while(abs(M)>abs(vthreshold))
        	event_location=[event_location;kloc iloc jc time];

        counter=counter+1;
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;
        temp=circshift(vjc,[kdelta idelta]);
        temp(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;
        vjc=circshift(temp,[-kdelta -idelta]);
        [M,I] = max(vjc(:));
        [kloc, iloc] = ind2sub(s,I);

        ufieldt	=circshift( ufieldt	,[kdelta idelta]);
        vfieldt	=circshift( vfieldt	,[kdelta idelta]);
        wfieldt	=circshift( wfieldt	,[kdelta idelta]);

    	dudxt=circshift( dudxt ,[kdelta idelta]);
        dvdxt=circshift( dvdxt ,[kdelta idelta]);
        dwdxt=circshift( dwdxt ,[kdelta idelta]);

        dudyt=circshift( dudyt ,[kdelta idelta]);
        dvdyt=circshift( dvdyt ,[kdelta idelta]);
        dwdyt=circshift( dwdyt ,[kdelta idelta]);

        dudzt=circshift( dudzt ,[kdelta idelta]);
        dvdzt=circshift( dvdzt ,[kdelta idelta]);
        dwdzt=circshift( dwdzt ,[kdelta idelta]);

        vozt	=circshift( vozt 	,[kdelta idelta]);
        woyt	=circshift( woyt 	,[kdelta idelta]);

	polyxt    =circshift(polyxt        ,[kdelta idelta]);
	polyyt    =circshift(polyyt        ,[kdelta idelta]);
	polyzt    =circshift(polyzt        ,[kdelta idelta]);

        un=	un	+ufieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        vn=	vn	-vfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        wn=	wn	+wfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudxn=dudxn        +dudxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdxn=dvdxn        -dvdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdxn=dwdxn        +dwdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        dudyn=dudyn        -dudyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdyn=dvdyn        +dvdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdyn=dwdyn        -dwdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        dudzn=dudzn        +dudzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdzn=dvdzn        -dvdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdzn=dwdzn        +dwdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        vozn=	vozn	+vozt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        woyn=	woyn	+woyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

	polyxn=  polyxn   +polyxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
	polyyn=  polyyn   -polyyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
	polyzn=  polyzn   +polyzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);


        ufieldt	=circshift( ufieldt,-[kdelta idelta]);
        vfieldt	=circshift( vfieldt,-[kdelta idelta]);
        wfieldt	=circshift( wfieldt,-[kdelta idelta]);

    	dudxt=circshift( dudxt,-[kdelta idelta]);
        dvdxt=circshift( dvdxt,-[kdelta idelta]);
        dwdxt=circshift( dwdxt,-[kdelta idelta]);

        dudyt=circshift( dudyt,-[kdelta idelta]);
        dvdyt=circshift( dvdyt,-[kdelta idelta]);
        dwdyt=circshift( dwdyt,-[kdelta idelta]);

        dudzt=circshift( dudzt,-[kdelta idelta]);
        dvdzt=circshift( dvdzt,-[kdelta idelta]);
        dwdzt=circshift( dwdzt,-[kdelta idelta]);

        vozt	=circshift( vozt 	,-[kdelta idelta]);
        woyt	=circshift( woyt 	,-[kdelta idelta]);

	polyxt    =circshift(polyxt        ,-[kdelta idelta]);
	polyyt    =circshift(polyyt        ,-[kdelta idelta]);
	polyzt    =circshift(polyzt        ,-[kdelta idelta]);
    end

end
counter
fc=sprintf("../data/conditionaln_jcond_%03d.mat",jcond);
mc=matfile(fc,'Writable',true);
mc.u=un./counter;
mc.v=vn./counter;
mc.w=wn./counter;

mc.dudx=dudxn./counter;
mc.dvdx=dvdxn./counter;
mc.dwdx=dwdxn./counter;

mc.dudy=dudyn./counter;
mc.dvdy=dvdyn./counter;
mc.dwdy=dwdyn./counter;

mc.dudz=dudzn./counter;
mc.dvdz=dvdzn./counter;
mc.dwdz=dwdzn./counter; 

mc.voz=vozn./counter;
mc.woy=woyn./counter;

mc.polyx=polyxn./counter;
mc.polyy=polyyn./counter;
mc.polyz=polyzn./counter;

