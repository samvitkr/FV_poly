jcond=188;
%mf=matfile('../data/lsevp_field_j_188.mat')

fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)

m1=matfile(fvgp,'Writable',true)
m2=matfile(fvgn,'Writable',true)

%fn=sprintf('../data/vel_corr_reflect_j_%03d.mat',jcond);
%mf=matfile(fn,"Writable",true)
%fn=sprintf('../data/lse_coeff_j_%03d.mat',jcond);
%mf=matfile(fn,"Writable",true)
%fn=sprintf('../data/mean_profiles.mat')
%mf=matfile(fn)
%m=matfile('../data/transferfields_0107000.mat')
%m=matfile('../data/velgradfield_lseQ4_j_156.mat')
%fpu=sprintf('../data/velgradfield_lseQ4_j_156.mat');
%mpu=matfile(fpu,'Writable',true)
%jcond=156;
%fnu=sprintf('../data/velgradfield_lseQ4_ufil_j_%03d.mat',jcond);
%mu=matfile(fnu,'Writable',true)
%fnu=sprintf('../data/velgradfield_lseQ4_dfil_j_%03d.mat',jcond);
%mu=matfile(fnu,'Writable',true)
%ft=sprintf('../data/velfield_lse_voz_j_188.mat')
%mt=matfile(ft)
%ft=sprintf('../data/velgradfx_voz_field_lseQ2_j_156.mat')
%m2=matfile(ft)
%ft=sprintf('../data/velgradfx_voz_field_lseQ4_j_156.mat')
%m4=matfile(ft)
%md=matfile('../data/velgrad_transfer_ddfilter_0040000.mat')
