function visualize_symmetry_lines(V,F,sampleIndex, sampleIndexSym)
%Hui Wang


%Plot the triangular mesh
Color0 = 0 * V + 211 / 255;%Grey color
trimesh(F,V(:,1),V(:,2),V(:,3), 'FaceVertexCData', Color0,'FaceColor','interp','EdgeColor', 'none','FaceAlpha',0.5);
set(gcf,'color','w')
axis off;
axis equal;
hold on;

%Plot the point-to-point correspondences of the symmetry
Xstart = V(sampleIndex,1);
Ystart = V(sampleIndex,2);
Zstart = V(sampleIndex,3);

Xend = V(sampleIndexSym,1);
Yend = V(sampleIndexSym,2);
Zend = V(sampleIndexSym,3);

plot3([Xstart, Xend]', [Ystart, Yend]', [Zstart, Zend]','Color',[0 0 1], 'LineWidth', 1.5);