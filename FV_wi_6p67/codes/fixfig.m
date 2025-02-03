%clim([-50 50])
clim([-0.1 0.1])
zlim([0.01 1])
%xlim([0 2*pi])
%xlim([3 6])
%ylim([ 3 6])
xlabel('x')
ylabel('z')
zlabel('y')
colormap redblue
view(-30,30)
camlight
lighting gouraud
%lightangle(-15,25)
lightangle(-15,-60)
OptionZ.FrameRate=100;OptionZ.Duration=20;OptionZ.Periodic=true;
%CaptureFigVid([-45,30;45,15;135,0;225,15;315,30], 'l2rms_flp_40000',OptionZ)
%CaptureFigVid([-45,30;45,15;135,0;225,15;315,30], 'l2rms_poly_unf',OptionZ)