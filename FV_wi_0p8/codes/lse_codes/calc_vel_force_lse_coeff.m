clear all
clc
Nx=512;
Nz=384;
Ny=220;
jcond=156;
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
fn=sprintf('vel_corr_j_%03d.mat',jcond);
%m=matfile(fn);
load(fn);

fm=sprintf('vel_force_corr_j_%03d.mat',jcond);
%m=matfile(fn);
load(fm);

jc=j;
uij=[Ruu(1,1,jc),Ruv(1,1,jc),Ruw(1,1,jc);...
    Rvu(1,1,jc),Rvv(1,1,jc),Rvw(1,1,jc);...
    Rwu(1,1,jc),Rwv(1,1,jc),Rww(1,1,jc)];
uij=uij.';

for j=1:Ny/2
    j
    for k=1:Nz
        for i=1:Nx
            Rij=[Ruu(k,i,j),Ruv(k,i,j),Ruw(k,i,j),Rufx(k,i,j);...
                Rvu(k,i,j),Rvv(k,i,j),Rvw(k,i,j),Rvfx(k,i,j);...
                Rwu(k,i,j),Rwv(k,i,j),Rww(k,i,j),Rwfx(k,i,j)];
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

        end
    end
end


fn=sprintf('vel_lse_j_%03d.mat',jcond);
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

mf.yCheb=yCheb;
mf.j=jc;



