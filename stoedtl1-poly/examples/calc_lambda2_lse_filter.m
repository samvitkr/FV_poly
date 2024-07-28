clear;  clc;  close all;
addpath '/home/skumar67/data-geyink1/skumar67/dmsuite'
sourceFolder='../';
% automatically add code modules to matlab path
functionCallStack = dbstack;
[scriptFolder, ~, ~] = fileparts(which(functionCallStack(1).file));
[sourceFolder, ~, ~] = fileparts(scriptFolder);
addpath(fullfile(sourceFolder, 'data_readers'), fullfile(sourceFolder, 'utilities'));

% user input: specify run parameters. This assumes your data is stored in a directory with full path dataFolder/runName
%runName = 'test_wi_6-67_1';
%fileNr = 79000;
ngx = 513;  % ng{i}: grid size of restart files in coordinate i, including periodic points in x, z.
ngy = 321;  %        ghost points do not have to be included here, they are accounted for in the data reading routines
ngz = 385;
nCheb = 220;  % number of Chebyshev collocation points in y (has to be determined a priori)
Lz = 2.0 * pi;
dataFolder = fullfile(getenv('HOME'), 'research_data', 'polymer');
Wi=0.8;
beta=0.9;
L_max=100;
a=1-3/L_max^2;
re=4667
% read grid and generate operators
%runFolder = fullfile(dataFolder, runName);
%runFolder= 'C:\Users\samvi\Dropbox\SimonsProject\Finite_Vol_data\Visco_data_sample';
%runFolder='/home/skumar67/data-geyink1/skumar67/FV_visco';
runFolder='/home/skumar67/data-geyink1/skumar67/FV_poly_codes/FV_wi_0p8/data'
[xGridPointsDns, yGridPointsDns, ~, ~] = read_grid(runFolder, ngx, ngy);
deltaX = xGridPointsDns(2) - xGridPointsDns(1);
deltaZ = Lz / (ngz - 1);
[yCheb, ~] = chebdif(nCheb, 1);
interpY = InterpolationOperatorY(yGridPointsDns, yCheb);
ft = SpatialFourierTransform();



jcond=156;
fn=sprintf('velfield_lse_uufilter_j_%03d.mat',jcond);
fn=fullfile(runFolder,fn);
mv=matfile(fn,'Writable',true)


Nx=512;
Ny=220;
Nz=384;
S_11	=single(zeros(Nz,Nx,Ny));
S_12	=single(zeros(Nz,Nx,Ny));
S_13	=single(zeros(Nz,Nx,Ny));
S_22	=single(zeros(Nz,Nx,Ny));
S_23	=single(zeros(Nz,Nx,Ny));
S_33	=single(zeros(Nz,Nx,Ny));
O_21	=single(zeros(Nz,Nx,Ny));
O_13	=single(zeros(Nz,Nx,Ny));
O_32	=single(zeros(Nz,Nx,Ny));

lambda2	=single(zeros(Nz,Nx,Ny));
Q	=single(zeros(Nz,Nx,Ny));


[uFourier, kx, kz] = ft.transform_to_fourier(mv.uQ2, deltaX, deltaZ);
vFourier = ft.transform_to_fourier(mv.vQ2);  % no need to pass sampling rate, wavenumber object is needed only once
wFourier = ft.transform_to_fourier(mv.wQ2);

uFourierq4 = ft.transform_to_fourier(mv.uQ4);
vFourierq4 = ft.transform_to_fourier(mv.vQ4);  % no need to pass sampling rate, wavenumber object is needed only once
wFourierq4 = ft.transform_to_fourier(mv.wQ4);



dx = FirstDerivativeXFourier(kx.get_radial_frequency_vector());
dz = FirstDerivativeZFourier(kz.get_radial_frequency_vector());
dy = FirstDerivativeYChebyshev(nCheb);

dudx=ft.transform_to_physical(dx.compute_derivative(uFourier));
dvdx=ft.transform_to_physical(dx.compute_derivative(vFourier));
dwdx=ft.transform_to_physical(dx.compute_derivative(wFourier));
dudy=ft.transform_to_physical(dy.compute_derivative(uFourier));
dvdy=ft.transform_to_physical(dy.compute_derivative(vFourier));
dwdy=ft.transform_to_physical(dy.compute_derivative(wFourier));
dudz=ft.transform_to_physical(dz.compute_derivative(uFourier));
dvdz=ft.transform_to_physical(dz.compute_derivative(vFourier));
dwdz=ft.transform_to_physical(dz.compute_derivative(wFourier));
	S_11=dudx;
        S_12=0.5*( dudy+dvdx );
        S_13=0.5*( dudz+dwdx );
        S_22=dvdy;
        S_23=0.5*( dwdy+dvdz );
        S_33=dwdz;
	omegaZ=dvdx-dudy;
	omegaY=dudz-dwdx;
	omegaX=dwdy-dvdz;
        O_21 = 0.5*omegaZ(:,:,:);
        O_13 = 0.5*omegaY(:,:,:);
        O_32 = 0.5*omegaX(:,:,:);
        O = zeros(3,3);
        S = zeros(3,3);
	for k =Ny/2+1:Ny
                for i =1:Nx
                for j =1:Nz
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
        %fl=sprintf("lambdaflp_%07d",time);
	%fl=sprintf("lambdafhp_%07d",time);
        %fl=sprintf("lambda_%07d",time);
	%ml=matfile(fl,'Writable',true);
        mv.Q2lambda2=single(lambda2);
        mv.Q2Q=single(Q);

dudx=ft.transform_to_physical(dx.compute_derivative(uFourierq4));
dvdx=ft.transform_to_physical(dx.compute_derivative(vFourierq4));
dwdx=ft.transform_to_physical(dx.compute_derivative(wFourierq4));
dudy=ft.transform_to_physical(dy.compute_derivative(uFourierq4));
dvdy=ft.transform_to_physical(dy.compute_derivative(vFourierq4));
dwdy=ft.transform_to_physical(dy.compute_derivative(wFourierq4));
dudz=ft.transform_to_physical(dz.compute_derivative(uFourierq4));
dvdz=ft.transform_to_physical(dz.compute_derivative(vFourierq4));
dwdz=ft.transform_to_physical(dz.compute_derivative(wFourierq4));
        S_11=dudx;
        S_12=0.5*( dudy+dvdx );
        S_13=0.5*( dudz+dwdx );
        S_22=dvdy;
        S_23=0.5*( dwdy+dvdz );
        S_33=dwdz;
        omegaZ=dvdx-dudy;
        omegaY=dudz-dwdx;
        omegaX=dwdy-dvdz;
        O_21 = 0.5*omegaZ(:,:,:);
        O_13 = 0.5*omegaY(:,:,:);
        O_32 = 0.5*omegaX(:,:,:);
        O = zeros(3,3);
        S = zeros(3,3);
	for k =Ny/2+1:Ny
                for i =1:Nx
                for j =1:Nz
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
	mv.Q4lambda2=single(lambda2);
        mv.Q4Q=single(Q);

