%% using the cmu.der package to compute derivatives
% John Kitchin

x = linspace(0,pi,100);
y = sin(x);

dydx = cmu.der.derc(y,x);

yp = cos(x); %analytical derivative

figure
hold on
plot(x,y,'k-')   %black line
plot(x,dydx)     % default blue line
plot(x,yp,'r-.') %red dot-dashed line
legend 'y' 'dydx' 'y'