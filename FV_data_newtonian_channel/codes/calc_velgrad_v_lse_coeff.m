clear all
clc
Nx=512;
Nz=384;
Ny=220;
jcond=194;
L11=zeros(Nz,Nx,Ny/2);
L12=zeros(Nz,Nx,Ny/2);
%L13=zeros(Nz,Nx,Ny/2);
L21=zeros(Nz,Nx,Ny/2);
L22=zeros(Nz,Nx,Ny/2);
%L23=zeros(Nz,Nx,Ny/2);
L31=zeros(Nz,Nx,Ny/2);
L32=zeros(Nz,Nx,Ny/2);
%L33=zeros(Nz,Nx,Ny/2);

L41=zeros(Nz,Nx,Ny/2);
L42=zeros(Nz,Nx,Ny/2);
L51=zeros(Nz,Nx,Ny/2);
L52=zeros(Nz,Nx,Ny/2);
L61=zeros(Nz,Nx,Ny/2);
L62=zeros(Nz,Nx,Ny/2);

L71=zeros(Nz,Nx,Ny/2);
L72=zeros(Nz,Nx,Ny/2);
L81=zeros(Nz,Nx,Ny/2);
L82=zeros(Nz,Nx,Ny/2);
L91=zeros(Nz,Nx,Ny/2);
L92=zeros(Nz,Nx,Ny/2);

L101=zeros(Nz,Nx,Ny/2);
L102=zeros(Nz,Nx,Ny/2);
L111=zeros(Nz,Nx,Ny/2);
L112=zeros(Nz,Nx,Ny/2);
L121=zeros(Nz,Nx,Ny/2);
L122=zeros(Nz,Nx,Ny/2);

L131=zeros(Nz,Nx,Ny/2);
L132=zeros(Nz,Nx,Ny/2);

fn=sprintf('../data/velgrad_corr_j_%03d.mat',jcond);

%m=matfile(fn);
load(fn);
jc=j;
uij=[Rvv(1,1,jc)];...,Rvw(1,1,jc)];...
    %Rwu(1,1,jc),Rwv(1,1,jc),Rww(1,1,jc)];
uij=uij.';

for j=1:Ny/2
    j
    for k=1:Nz
        for i=1:Nx
            Rij=[Rvu(k,i,j),Rvv(k,i,j),Rvw(k,i,j),...
		Rvdudx(k,i,j),Rvdvdx(k,i,j),Rvdwdx(k,i,j),...
                Rvdudy(k,i,j),Rvdvdy(k,i,j),Rvdwdy(k,i,j),...
                Rvdudz(k,i,j),Rvdvdz(k,i,j),Rvdwdz(k,i,j),Rvfx(k,i,j)];
               % Rwu(k,i,j),Rwv(k,i,j),Rww(k,i,j)];
            Rij=Rij.';
            L=Rij*(inv(uij));
            L11(k,i,j)=L(1,1);
            L21(k,i,j)=L(2,1);
            L31(k,i,j)=L(3,1);
            L41(k,i,j)=L(4,1);
	    L51(k,i,j)=L(5,1);
	    L61(k,i,j)=L(6,1);
	    L71(k,i,j)=L(7,1);
	    L81(k,i,j)=L(8,1);
	    L91(k,i,j)=L(9,1);
	    L101(k,i,j)=L(10,1);
            L111(k,i,j)=L(11,1);
            L121(k,i,j)=L(12,1);
	    L131(k,i,j)=L(13,1);

        end
    end
end


fn=sprintf('../data/velgrad_v_lse_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.L11=real(L11);
mf.L21=real(L21);
mf.L31=real(L31);
mf.yCheb=yCheb;
mf.j=jc;
mf.L41=real(L41);
mf.L51=real(L51);
mf.L61=real(L61);
mf.L71=real(L71);
mf.L81=real(L81);
mf.L91=real(L91);
mf.L101=real(L101);
mf.L111=real(L111);
mf.L121=real(L121);
mf.L131=real(L131);
