clear
clc
load('../data/ygrid.mat')
y=yCheb+1;
m2=matfile('../data/velgradfield_lseQ2_voz_j_156.mat');
m4=matfile('../data/velgradfield_lseQ4_voz_j_156.mat');
nl2=squeeze(mean( m2.v.*(m2.dvdx-m2.dudy) - m2.w.*(m2.dudz-m2.dwdx) ,[1 2])) ; 
nl4=squeeze(mean( m4.v.*(m4.dvdx-m4.dudy) - m4.w.*(m4.dudz-m4.dwdx) ,[1 2])) ; 

voz2=squeeze(mean( m2.v.*(m2.dvdx-m2.dudy)  ,[1 2])) ; 
voz4=squeeze(mean( m4.v.*(m4.dvdx-m4.dudy) ,[1 2])) ; 

woy2=squeeze(mean(  m2.w.*(m2.dudz-m2.dwdx) ,[1 2])) ; 
woy4=squeeze(mean(  m4.w.*(m4.dudz-m4.dwdx) ,[1 2])) ; 
%%
y=yCheb(111:end)+1;

hold on
plot(-nl2,y,'-r')
plot(-nl4,y,'-b')
plot(-voz2,y,'--r')
plot(-voz4,y,'--b')
plot(woy2,y,':r','LineWidth',1.5)
plot(woy4,y,':b','LineWidth',1.5)
hold off
yline(y(46))
grid on
xlabel('(v\omega_z-w\omega_y)/(-U_b^3/H)')
ylabel('y')
legend('v<0','v>0')