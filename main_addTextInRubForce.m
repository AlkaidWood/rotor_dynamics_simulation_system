ymax = 6000;
ymin = 0;
x_point = [6e-3 9.8e-3  12e-3];
y_point = ymax;
xline = [x_point(1),x_point(1),x_point(2),x_point(2);
         x_point(2),x_point(2),x_point(3),x_point(3)];
y1ine = [ ymin, y_point, y_point, ymin;
         -ymin,-y_point,-y_point,-ymin;];

blockcolor1=[254, 239, 237]/255;
blockcolor2=[240, 244, 246]/255;
blockcolor3=[236, 225, 212]/255;

fill_1 = fill(xline(1,:),[0 ymax ymax 0],blockcolor2, 'HandleVisibility','off');hold on
fill_2 = fill(xline(2,:),[0 ymax ymax 0],blockcolor1, 'HandleVisibility','off');

set(fill_1,'edgecolor','none');
set(fill_2,'edgecolor','none');
uistack(fill_1, 'bottom')
uistack(fill_2, 'bottom')

set(gca, 'layer', 'top');

textHeight = 5000;
text(7.9e-3, textHeight,'无碰摩区域','VerticalAlignment','bottom', 'HorizontalAlignment','center','Fontname', '宋体','FontSize',8);
text(11e-3, textHeight,'碰摩区域','VerticalAlignment','bottom', 'HorizontalAlignment','center','Fontname', '宋体','FontSize',8);