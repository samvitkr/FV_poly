clear all
clc
clear
%numWorkers = 16; % Use 48 processors on the current node

% Start the parallel pool with the specified number of workers
%parpool('local', numWorkers);
Nx=512;
Nz=384;
Ny=220;
jcond=156;
%jc=Ny-jcond+1;

fn=sprintf('../data/lse_coeff_B_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);

L11=		single(zeros(Nz,Nx,Ny/2));
L21=		single(zeros(Nz,Nx,Ny/2));
L31=		single(zeros(Nz,Nx,Ny/2));
L41=		single(zeros(Nz,Nx,Ny/2));
L51=		single(zeros(Nz,Nx,Ny/2));
L61=		single(zeros(Nz,Nx,Ny/2));
L71=		single(zeros(Nz,Nx,Ny/2));
L81=		single(zeros(Nz,Nx,Ny/2));
L91=		single(zeros(Nz,Nx,Ny/2));
L101=		single(zeros(Nz,Nx,Ny/2));
L111=		single(zeros(Nz,Nx,Ny/2));
L121=		single(zeros(Nz,Nx,Ny/2));
L131=		single(zeros(Nz,Nx,Ny/2));
L141=		single(zeros(Nz,Nx,Ny/2));
L151=           single(zeros(Nz,Nx,Ny/2));
L161=           single(zeros(Nz,Nx,Ny/2));
L171=           single(zeros(Nz,Nx,Ny/2));

fn=sprintf('../data/vel_corr_reflect_j_%03d.mat',jcond);
m=matfile(fn);
uij=[m.Rvv(1,1,jcond-Ny/2)];
uij=uij.';
clear m

fn=sprintf('../data/vel_corr_reflect_B_j_%03d.mat',jcond);
load(fn)

for j=1:Ny/2
    j
    	for k=1:Nz
        for i=1:Nx
            Rij=[Rvu(k,i,j),Rvv(k,i,j),Rvw(k,i,j),...
		Rvdudx(k,i,j),Rvdvdx(k,i,j),Rvdwdx(k,i,j),...
                Rvdudy(k,i,j),Rvdvdy(k,i,j),Rvdwdy(k,i,j),...
                Rvdudz(k,i,j),Rvdvdz(k,i,j),Rvdwdz(k,i,j),...
		Rvvoz(k,i,j),Rvwoy(k,i,j),...
		Rvfx(k,i,j),Rvfy(k,i,j),Rvfz(k,i,j)];
            Rij=Rij.';
            L=Rij/(uij);
            L11(k,i,j)=	single(real(L(1,1)));
            L21(k,i,j)=	single(real(L(2,1)));
            L31(k,i,j)=	single(real(L(3,1)));
            L41(k,i,j)=	single(real(L(4,1)));
	    L51(k,i,j)=	single(real(L(5,1)));
	    L61(k,i,j)=	single(real(L(6,1)));
	    L71(k,i,j)=	single(real(L(7,1)));
	    L81(k,i,j)=	single(real(L(8,1)));
	    L91(k,i,j)=	single(real(L(9,1)));
	    L101(k,i,j)=	single(real(L(10,1)));
            L111(k,i,j)=	single(real(L(11,1)));
            L121(k,i,j)=	single(real(L(12,1)));
	    L131(k,i,j)=	single(real(L(13,1)));
	    L141(k,i,j)=	single(real(L(14,1)));
	    L151(k,i,j)=  single(real(L(15,1)));
            L161(k,i,j)=  single(real(L(16,1)));
            L171(k,i,j)=  single(real(L(17,1)));

        end
    	end
end
	 mf.L11=L11; 
         mf.L21=L21; 
         mf.L31=L31; 
         mf.L41=L41; 
         mf.L51=L51; 
         mf.L61=L61; 
         mf.L71=L71; 
         mf.L81=L81; 
         mf.L91=L91; 
        mf.L101=L101; 
        mf.L111=L111; 
        mf.L121=L121; 
        mf.L131=L131; 
        mf.L141=L141;
	mf.L151=L151;
        mf.L161=L161;
        mf.L171=L171;
