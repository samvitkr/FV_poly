
x = 0:.1:1;
A = exp(x);
x1=180;
y1=180;
x2=400;
y2=500;

figure('OuterPosition',...
    [x1 y1 x2 y2]);
plot(x,A);
print('outfile','-dpng');

