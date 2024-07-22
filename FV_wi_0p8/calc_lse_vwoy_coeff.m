    clear all
clc
Nx=512;
Nz=384;
Ny=220;
jcond=188;

L11=zeros(Nz,Nx,Ny/2);
L12=zeros(Nz,Nx,Ny/2);
L13=zeros(Nz,Nx,Ny/2);

L21=zeros(Nz,Nx,Ny/2);
L22=zeros(Nz,Nx,Ny/2);
L23=zeros(Nz,Nx,Ny/2);

L31=zeros(Nz,Nx,Ny/2);
L32=zeros(Nz,Nx,Ny/2);
L33=zeros(Nz,Nx,Ny/2);

L41=zeros(Nz,Nx,Ny/2);
L42=zeros(Nz,Nx,Ny/2);
L43=zeros(Nz,Nx,Ny/2);

L51=zeros(Nz,Nx,Ny/2);
L52=zeros(Nz,Nx,Ny/2);
L53=zeros(Nz,Nx,Ny/2);

L61=zeros(Nz,Nx,Ny/2);
L62=zeros(Nz,Nx,Ny/2);
L63=zeros(Nz,Nx,Ny/2);

%mn=matfile('transferfields_mean.mat');
%fn=sprintf('voz_woy_corr_j_%03d.mat',jcond);
%m=matfile(fn);
fn=sprintf('vwoy_corr_j_%03d.mat',jcond);
load(fn);

%fm=sprintf('vel_force_corr_j_%03d.mat',jcond);
%m=matfile(fn);
%load(fm);

jc=j;
%woy=mean(mn.woy(:,:,jc),'all');
%uij=[	Rww(1,1,jc),woy ,Rwv(1,1,jc);...
%    	woy	   ,oyoy,voy;...
%	Rvw(1,1,jc),voy ,Rvv(1,1,jc)];
uij=[vv  vw  voy;...
     vw  ww  woy;...
     voy woy oyoy]; 
uij=uij.';

for j=1:Ny/2
    j
    for k=1:Nz
        for i=1:Nx
            Rij=[ Rvu(k,i,j), Rvv(k,i,j), Rvw(k,i,j), Rvfx(k,i,j), Rvfy(k,i,j), Rvfz(k,i,j); 
	    	  Rwu(k,i,j), Rwv(k,i,j), Rww(k,i,j), Rwfx(k,i,j), Rwfy(k,i,j), Rwfz(k,i,j);...
                 Royu(k,i,j),Royv(k,i,j),Royw(k,i,j),Royfx(k,i,j),Royfy(k,i,j),Royfz(k,i,j)];
	  
            Rij=Rij.';
       
       	    L=Rij*(inv(uij));

            L11(k,i,j)=L(1,1);
            L12(k,i,j)=L(1,2);
            L13(k,i,j)=L(1,3);

            L21(k,i,j)=L(2,1);
            L22(k,i,j)=L(2,2);
            L23(k,i,j)=L(2,3);

            L31(k,i,j)=L(3,1);
            L32(k,i,j)=L(3,2);
            L33(k,i,j)=L(3,3);

	    L41(k,i,j)=L(4,1);
            L42(k,i,j)=L(4,2);
            L43(k,i,j)=L(4,3);

	    L51(k,i,j)=L(5,1);
            L52(k,i,j)=L(5,2);
            L53(k,i,j)=L(5,3);

	    L61(k,i,j)=L(6,1);
            L62(k,i,j)=L(6,2);
            L63(k,i,j)=L(6,3);
        end
    end
end


fn=sprintf('vwoy_lse_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);

mf.L11=real(L11);
mf.L12=real(L12);
mf.L13=real(L13);

mf.L21=real(L21);
mf.L22=real(L22);
mf.L23=real(L23);

mf.L31=real(L31);
mf.L32=real(L32);
mf.L33=real(L33);

mf.L41=real(L41);
mf.L42=real(L42);
mf.L43=real(L43);

mf.L51=real(L51);
mf.L52=real(L52);
mf.L53=real(L53);

mf.L61=real(L61);
mf.L62=real(L62);
mf.L63=real(L63);

mf.j=jc;
