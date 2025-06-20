jcond=171;
%mf=matfile('../data/lsevp_field_j_188.mat')

fvgp=sprintf('../data/conditionalp_jcond_%03d.mat',jcond)
fvgn=sprintf('../data/conditionaln_jcond_%03d.mat',jcond)

mp=matfile(fvgp,'Writable',true)
mn=matfile(fvgn,'Writable',true)

