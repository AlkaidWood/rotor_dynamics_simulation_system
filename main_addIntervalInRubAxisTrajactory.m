interval = 6.55e-3;
r = interval;
a = 0;
b = 0;
theta = 0:pi/50:2*pi;
x = a + r*cos(theta);
y = b + r*sin(theta);

mycolor = [ 0.40784,0.5804,0.651;
            0.94902,0.37647,0.32157;
            191/255,155/255,111/255;
            242/255,203/255,5/255];
        
        
h = gcf;
figure(h);hold on
plot(x, y, '-','LineWidth',1,'color',mycolor(2,:));

