Nx=512;
Ny=220;
Nz=384;
jcond=156;
S_11	=single(zeros(Nz,Nx,Ny/2));
S_12	=single(zeros(Nz,Nx,Ny/2));
S_13	=single(zeros(Nz,Nx,Ny/2));
S_22	=single(zeros(Nz,Nx,Ny/2));
S_23	=single(zeros(Nz,Nx,Ny/2));
S_33	=single(zeros(Nz,Nx,Ny/2));
O_21	=single(zeros(Nz,Nx,Ny/2));
O_13	=single(zeros(Nz,Nx,Ny/2));
O_32	=single(zeros(Nz,Nx,Ny/2));
lambda2	=single(zeros(Nz,Nx,Ny/2));
Q	=single(zeros(Nz,Nx,Ny/2));
fvg=sprintf('../data/velgrad_voz_field_lseQ4_j_%03d.mat',jcond);
mvg=matfile(fvg,'Writable',true)
fl=size(mvg.dudy);
if(fl(3)>Ny/2)
	mvg.u=mvg.u(1:end,1:end,111:end);
	mvg.v=mvg.v(1:end,1:end,111:end);
	mvg.w=mvg.w(1:end,1:end,111:end);
	mvg.dudx=mvg.dudx(1:end,1:end,111:end);
	mvg.dvdx=mvg.dvdx(1:end,1:end,111:end);
	mvg.dwdx=mvg.dwdx(1:end,1:end,111:end);
	mvg.dudy=mvg.dudy(1:end,1:end,111:end);
        mvg.dvdy=mvg.dvdy(1:end,1:end,111:end);
        mvg.dwdy=mvg.dwdy(1:end,1:end,111:end);
	mvg.dudz=mvg.dudz(1:end,1:end,111:end);
        mvg.dvdz=mvg.dvdz(1:end,1:end,111:end);
        mvg.dwdz=mvg.dwdz(1:end,1:end,111:end);

end
S_11=mvg.dudx;
S_12=0.5*( mvg.dudy+mvg.dvdx );
S_13=0.5*( mvg.dudz+mvg.dwdx );
S_22=mvg.dvdy;
S_23=0.5*( mvg.dwdy+mvg.dvdz );
S_33=mvg.dwdz;

omegaZ=mvg.dvdx-mvg.dudy;
omegaY=mvg.dudz-mvg.dwdx;
omegaX=mvg.dwdy-mvg.dvdz;
O_21 = 0.5*omegaZ(:,:,:);
O_13 = 0.5*omegaY(:,:,:);
O_32 = 0.5*omegaX(:,:,:);

O = zeros(3,3);
S = zeros(3,3);

for j =1:Nz
        for i =1:Nx
        for k =1:Ny/2
        S(1,1) = S_11(j,i,k);%mvg.dudx(i,j,k);
        S(1,2) = S_12(j,i,k);%0.5*( mvg.dudy(i,j,k) +mvg.dvdx(i,j,k) );
        S(1,3) = S_13(j,i,k);%0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
        S(2,1) = S(1,2);
        S(2,2) = S_22(j,i,k);% mvg.dvdy(i,j,k);
        S(2,3) = S_23(j,i,k);%0.5*( mvelgz.dvdz(i,j,kstart+k)+mvg.dwdy(i,j,k));
        S(3,1) = S(1,3);%,ks 0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
        S(3,2) = S(2,3);
        S(3,3) = S_33(j,i,k);%mvelgz.dwdz(i,j,kstart+k);

        O(1,3) = O_13(j,i,k);%0.5*mo.omega_y(i,j,kstart+k);
        O(2,1) = O_21(j,i,k);%0.5*mo.omega_z(i,j,kstart+k);
        O(3,2) = O_32(j,i,k);% 0.5*mo.omega_x(i,j,kstart+k);
        O(1,2) =-O(2,1);
        O(2,3) =-O(3,2);
        O(3,1) =-O(1,3);

        A = S*S + O*O;
        B = O*O';
	C = S*S';

        ll = sort(eig(A));
        lambda2(j,i,k) = single(ll(2));
        Q(j,i,k)= 0.5*(trace(B)-trace(C));
        end
        end
end
mvg.lambda2=single(lambda2);
mvg.Q=single(Q);

%%
%clear mvg S_11 S_12 S_13 S_22 S_23 S_33 O_21 O_13 O_32
%fvg=sprintf('../data/velgradfx_voz_field_lseQ4_j_%03d.mat',jcond)
%mvg=matfile(fvg,'Writable',true);
%fl=size(mvg.dudy);
%if(fl(3)>Ny/2)
%        mvg.u=mvg.u(1:end,1:end,111:end);
%        mvg.v=mvg.v(1:end,1:end,111:end);
%        mvg.w=mvg.w(1:end,1:end,111:end);
%        mvg.dudx=mvg.dudx(1:end,1:end,111:end);
%        mvg.dvdx=mvg.dvdx(1:end,1:end,111:end);
%        mvg.dwdx=mvg.dwdx(1:end,1:end,111:end);
%        mvg.dudy=mvg.dudy(1:end,1:end,111:end);
%        mvg.dvdy=mvg.dvdy(1:end,1:end,111:end);
%        mvg.dwdy=mvg.dwdy(1:end,1:end,111:end);
%        mvg.dudz=mvg.dudz(1:end,1:end,111:end);
%        mvg.dvdz=mvg.dvdz(1:end,1:end,111:end);
%        mvg.dwdz=mvg.dwdz(1:end,1:end,111:end);
%
%end
%S_11=mvg.dudx;
%S_12=0.5*( mvg.dudy+mvg.dvdx );
%S_13=0.5*( mvg.dudz+mvg.dwdx );
%S_22=mvg.dvdy;
%S_23=0.5*( mvg.dwdy+mvg.dvdz );
%S_33=mvg.dwdz;
%
%omegaZ=mvg.dvdx-mvg.dudy;
%omegaY=mvg.dudz-mvg.dwdx;
%omegaX=mvg.dwdy-mvg.dvdz;
%O_21 = 0.5*omegaZ(:,:,:);
%O_13 = 0.5*omegaY(:,:,:);
%O_32 = 0.5*omegaX(:,:,:);
%
%O = zeros(3,3);
%S = zeros(3,3);
%
%for j =1:Nz
%        for i =1:Nx
%        for k =1:Ny/2
%        S(1,1) = S_11(j,i,k);%mvg.dudx(i,j,k);
%        S(1,2) = S_12(j,i,k);%0.5*( mvg.dudy(i,j,k) +mvg.dvdx(i,j,k) );
%        S(1,3) = S_13(j,i,k);%0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
%        S(2,1) = S(1,2);
%        S(2,2) = S_22(j,i,k);% mvg.dvdy(i,j,k);
%        S(2,3) = S_23(j,i,k);%0.5*( mvelgz.dvdz(i,j,kstart+k)+mvg.dwdy(i,j,k));
%        S(3,1) = S(1,3);%,ks 0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
%        S(3,2) = S(2,3);
%        S(3,3) = S_33(j,i,k);%mvelgz.dwdz(i,j,kstart+k);
%
%        O(1,3) = O_13(j,i,k);%0.5*mo.omega_y(i,j,kstart+k);
%        O(2,1) = O_21(j,i,k);%0.5*mo.omega_z(i,j,kstart+k);
%        O(3,2) = O_32(j,i,k);% 0.5*mo.omega_x(i,j,kstart+k);
%        O(1,2) =-O(2,1);
%        O(2,3) =-O(3,2);
%        O(3,1) =-O(1,3);
%
%        A = S*S + O*O;
%        B = O*O';
%        C = S*S';
%
%        ll = sort(eig(A));
%        lambda2(j,i,k) = single(ll(2));
%        Q(j,i,k)= 0.5*(trace(B)-trace(C));
%        end
%        end
%end
%mvg.lambda2=single(lambda2);
%mvg.Q=single(Q);





