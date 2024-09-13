m=matfile("../data/velfield_lse_voz_j_156.mat",'Writable',true)
uQ2=m.uQ2+flip(m.uQ2,3);
wQ2=m.wQ2+flip(m.wQ2,3);
vQ2=m.vQ2-flip(m.vQ2,3);
fxQ2=m.fxQ2+flip(m.fxQ2,3);

uQ4=m.uQ4+flip(m.uQ4,3);
wQ4=m.wQ4+flip(m.wQ4,3);
vQ4=m.vQ4-flip(m.vQ4,3);
fxQ4=m.fxQ4+flip(m.fxQ4,3);

if( (sum(abs(m.uQ2(:,:,2)),'all')<1e-16)) 
    m.uQ2=uQ2;
m.vQ2=vQ2;
m.wQ2=wQ2;
end
if( (sum(abs(m.uQ4(:,:,2)),'all')<1e-16)) 
m.uQ4=uQ4;
m.vQ4=vQ4;
m.wQ4=wQ4;
end