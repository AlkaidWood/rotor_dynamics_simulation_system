% % 创建 textarrow
% annotation('textarrow',[0.222873900293255 0.274581500488759],...
% [0.69620253164557 0.707181188808324],...
% 'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
% 'String','$f_L$',...
% 'Interpreter','latex',...
% 'HeadWidth',3,...
% 'HeadLength',3,...
% 'FontSize',7);

% % 创建 textarrow
% annotation('textarrow',[0.777126099706745 0.747800586510264],...
% [0.405063291139241 0.261603375527426],...
% 'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
% 'String','$\frac{3}{2}f_L$',...
% 'Interpreter','latex',...
% 'HeadWidth',5,...
% 'HeadLength',3,...
% 'FontSize',7);

% % 创建 textarrow
% annotation('textarrow',[0.721407624633431 0.724340175953079],...
%     [0.440928270042194 0.293248945147679],...
%     'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
%     'String','$2f_L$',...
%     'Interpreter','latex',...
%     'HeadWidth',5,...
%     'HeadLength',3,...
%     'FontSize',7);
% 
% % 创建 textarrow
% annotation('textarrow',[0.782258064516129 0.783724340175955],...
%     [0.413502109704641 0.337552742616033],...
%     'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
%     'String','$3f_L$',...
%     'Interpreter','latex',...
%     'HeadWidth',5,...
%     'HeadLength',3,...
%     'FontSize',7);

% % 创建 textarrow
% annotation('textarrow',[0.175219941348974 0.226927541544478],...
%     [0.69620253164557 0.707181188808324],...
%     'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
%     'String','$f_L$',...
%     'Interpreter','latex',...
%     'HeadWidth',3,...
%     'HeadLength',3,...
%     'FontSize',7);

% % 创建 textarrow
% annotation('textarrow',[0.279166666666667 0.278980327468232],...
%     [0.560949298813377 0.471949121297774],...
%     'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
%     'String','$f_H$',...
%     'Interpreter','latex',...
%     'HeadWidth',3,...
%     'HeadLength',3,...
%     'FontSize',7);



% 创建 松动故障典型频率区域并标注
blockcolor1=[254, 239, 237]/255;
blockcolor2=[242, 96, 82]/255;
ymax = 10e7;
ymin = 0;
x_point = [30 100];
y_point = ymax;
xline = [x_point(1),x_point(1),x_point(2),x_point(2)];
y1ine = [ ymin, y_point, y_point, ymin;
         -ymin,-y_point,-y_point,-ymin;];
fill_1 = fill(xline(1,:),[0 ymax ymax 0],blockcolor2, 'HandleVisibility','off');
set(fill_1,'edgecolor','none');
uistack(fill_1, 'bottom')

axis = get(gca);
% yPos = 10e7;
yPos = 1500;
zPos = abs((axis.ZLim(2) - axis.ZLim(1)));
xPos = (x_point(2)+x_point(1))/2;
text(xPos, yPos, zPos*0.05,'$f_b$','rotation', 30,'VerticalAlignment','bottom', 'HorizontalAlignment','center','FontSize',7,'interpreter','latex');











