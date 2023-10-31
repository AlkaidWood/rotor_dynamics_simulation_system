

ymax = 4e-3;
ymin = 2.17e-3;
x_point = [0 2.8 4 6];
y_point = ymax;
xline = [x_point(1),x_point(1),x_point(2),x_point(2);
         x_point(2),x_point(2),x_point(3),x_point(3);
         x_point(3),x_point(3),x_point(4),x_point(4)];
y1ine = [ ymin, y_point, y_point, ymin;
         -ymin,-y_point,-y_point,-ymin;];

blockcolor1=[254, 239, 237]/255;
blockcolor2=[240, 244, 246]/255;
blockcolor3=[236, 225, 212]/255;
h = gcf;
figure(h); hold on

fill_1 = fill(xline(1,:),[0 ymax ymax 0],blockcolor2, 'HandleVisibility','off');hold on
fill_2 = fill(xline(1,:),[0 -ymax -ymax 0],blockcolor2, 'HandleVisibility','off');hold on
fill_3 = fill(xline(2,:),y1ine(1,:),blockcolor1, 'HandleVisibility','off');hold on
fill_4 = fill(xline(2,:),y1ine(2,:),blockcolor1, 'HandleVisibility','off');hold on
fill_5 = fill(xline(3,:),y1ine(1,:),blockcolor3, 'HandleVisibility','off');hold on
fill_6 = fill(xline(3,:),y1ine(2,:),blockcolor3, 'HandleVisibility','off');

set(fill_1,'edgecolor','none');
set(fill_2,'edgecolor','none');
set(fill_3,'edgecolor','none');
set(fill_4,'edgecolor','none');
set(fill_5,'edgecolor','none');
set(fill_6,'edgecolor','none');
uistack(fill_1, 'bottom')
uistack(fill_2, 'bottom')
uistack(fill_3, 'bottom')
uistack(fill_4, 'bottom')
uistack(fill_5, 'bottom')
uistack(fill_6, 'bottom')

set(gca, 'layer', 'top');
set(gcf,'Units','centimeters','Position',[6 6 14.8 4]);%Set the size of figure(for A4)

textHeight = 2.9e-3;
text(1.4, textHeight,'瞬态加速段（无碰摩）','VerticalAlignment','bottom', 'HorizontalAlignment','center','Fontname', '宋体','FontSize',8);
text(1.4, -textHeight,'瞬态加速段（无碰摩）','VerticalAlignment','top', 'HorizontalAlignment','center', 'Fontname', '宋体','FontSize',8);
text(3.4, textHeight,'全周碰摩','VerticalAlignment','bottom', 'HorizontalAlignment','center', 'Fontname', '宋体','FontSize',8);
text(3.4, -textHeight,'全周碰摩','VerticalAlignment','top', 'HorizontalAlignment','center', 'Fontname', '宋体','FontSize',8);
text(5, textHeight,'局部碰摩','VerticalAlignment','bottom', 'HorizontalAlignment','center', 'Fontname', '宋体','FontSize',8);
text(5, -textHeight,'局部碰摩','VerticalAlignment','top', 'HorizontalAlignment','center', 'Fontname', '宋体','FontSize',8);
% print(h, ['G:/大学硕士/毕业论文/论文/template_master-master/figures/rub/bifurcation/', 'Node-8-DOF-1'], '-dpng','-r400')