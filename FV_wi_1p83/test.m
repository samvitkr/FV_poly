m=matfile("velgrad_transfer_flp_0070000.mat")
muf=matfile("velgrad_0070000.mat")
%mt=matfile('testdudx.mat','Writable',true)
%mt.dudx_flp=m.dudx;
%mt.dudx=muf.dudx;
size(muf.dudxF)
