clear all
Nx=512;
Nz=384;
Ny=220;
jcond=188;
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
%L43=zeros(Nz,Nx,Ny/2);

%mn=matfile('transferfields_mean.mat');
fn=sprintf('../data/voz_corr_ddfilter_j_%03d.mat',jcond);
mdd=matfile(fn)
load(fn);

%fm=sprintf('vel_force_corr_j_%03d.mat',jcond);
%m=matfile(fn);
%load(fm);

jc=j;
voz=Rozv(1,1,jc);%mean(mn.voz(:,:,jc),'all');
uij=[Rvv(1,1,jc),	voz;...
     Rozv(1,1,jc),      ozoz];
uij=uij.';

for j=1:Ny/2
    j
    for k=1:Nz
        for i=1:Nx
	Rij=[Rvu(k,i,j),Rvv(k,i,j),Rvw(k,i,j),Rvfx(k,i,j);...
                Rozu(k,i,j),Rozv(k,i,j),Rozw(k,i,j),Rozfx(k,i,j)];
            Rij=Rij.';
            L=Rij*(inv(uij));
            L11(k,i,j)=L(1,1);
            L12(k,i,j)=L(1,2);
            %L13(k,i,j)=L(1,3);

            L21(k,i,j)=L(2,1);
            L22(k,i,j)=L(2,2);
            %L23(k,i,j)=L(2,3);

            L31(k,i,j)=L(3,1);
            L32(k,i,j)=L(3,2);
            %L33(k,i,j)=L(3,3);

	    L41(k,i,j)=L(4,1);
            L42(k,i,j)=L(4,2);
            %L43(k,i,j)=L(4,3);

        end
    end
end
fn=sprintf('../data/voz_lse_ddfilter_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.L11=real(L11);
mf.L12=real(L12);
%mf.L13=real(L13);
mf.L21=real(L21);
mf.L22=real(L22);
%mf.L23=real(L23);
mf.L31=real(L31);
mf.L32=real(L32);
%mf.L33=real(L33);
mf.L41=real(L41);
mf.L42=real(L42);
%mf.L43=real(L43);
mf.j=jc;


clear all
Nx=512;
Nz=384;
Ny=220;
jcond=188;
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
%L43=zeros(Nz,Nx,Ny/2);

%mn=matfile('transferfields_mean.mat');
fn=sprintf('../data/voz_corr_uufilter_j_%03d.mat',jcond);
muu=matfile(fn)
load(fn);

%fm=sprintf('vel_force_corr_j_%03d.mat',jcond);
%m=matfile(fn);
%load(fm);

jc=j;
voz=Rozv(1,1,jc);%mean(mn.voz(:,:,jc),'all');
uij=[Rvv(1,1,jc),       voz;...
     Rozv(1,1,jc),      ozoz];
uij=uij.';

for j=1:Ny/2
    j
    for k=1:Nz
        for i=1:Nx
		Rij=[Rvu(k,i,j),Rvv(k,i,j),Rvw(k,i,j),Rvfx(k,i,j);...
                Rozu(k,i,j),Rozv(k,i,j),Rozw(k,i,j),Rozfx(k,i,j)];
            Rij=Rij.';
            L=Rij*(inv(uij));
            L11(k,i,j)=L(1,1);
            L12(k,i,j)=L(1,2);
            %L13(k,i,j)=L(1,3);

            L21(k,i,j)=L(2,1);
            L22(k,i,j)=L(2,2);
            %L23(k,i,j)=L(2,3);

            L31(k,i,j)=L(3,1);
            L32(k,i,j)=L(3,2);
            %L33(k,i,j)=L(3,3);

            L41(k,i,j)=L(4,1);
            L42(k,i,j)=L(4,2);
            %L43(k,i,j)=L(4,3);

        end
    end
end
fn=sprintf('../data/voz_lse_uufilter_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.L11=real(L11);
mf.L12=real(L12);
%mf.L13=real(L13);
mf.L21=real(L21);
mf.L22=real(L22);
%mf.L23=real(L23);
mf.L31=real(L31);
mf.L32=real(L32);
%mf.L33=real(L33);
mf.L41=real(L41);
mf.L42=real(L42);
%mf.L43=real(L43);
mf.j=jc;

