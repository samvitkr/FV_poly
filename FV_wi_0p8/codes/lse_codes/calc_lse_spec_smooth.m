
win=3;
jcond=156;
fpu=sprintf('../data/velgrad_lse_spec_j_%03d.mat',jcond);
mpu=matfile(fpu,'Writable',true)
phiconv=mpu.phinl2;
[Nz Nx Ny] = size(phiconv);
phismooth=phiconv.*0;
for j =1:Ny
	val=spec_smoother(phiconv(:,:,j),win);
	phismooth(:,:,j)=val(1:Nz,1:Nx);
end
mpu.phinl2smooth=phismooth;

phiconv=mpu.phinl4;
[Nz Nx Ny] = size(phiconv);
phismooth=phiconv.*0;
for j =1:Ny
        val=spec_smoother(phiconv(:,:,j),win);
        phismooth(:,:,j)=val(1:Nz,1:Nx);
end
mpu.phinl4smooth=phismooth;
