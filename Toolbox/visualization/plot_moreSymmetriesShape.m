function plot_moreSymmetriesShape(V, F, originalColor, symmetryOrbit, CC, symmetryIndex)
%Hui Wang

figure;
set(gcf,'color','w','name','The original color and symmetry orbit')
trimesh(F,V(:,1),V(:,2),V(:,3), 'FaceVertexCData', originalColor, 'FaceColor', 'interp', 'EdgeColor', 'none')
view([-2 42])
axis off
axis equal;
hold on
plot3(V(symmetryOrbit,1), V(symmetryOrbit,2), V(symmetryOrbit,3), 'r*')

a = get(0);
figure('position',a.MonitorPositions);
set(gcf,'color','w','name','The computed symmetries')
numberOfSymmetry = size(symmetryIndex,2);
for i = 1:numberOfSymmetry
    subplot(2, 4, i);
    trimesh(F,V(:,1),V(:,2),V(:,3), 'FaceVertexCData', originalColor(symmetryIndex(:,i),:), 'FaceColor', 'interp', 'EdgeColor', 'none')
    view([-2 42])
    axis off
    axis equal;
    title(strcat('g_', num2str(i)))
end

a = get(0);
figure('position',a.MonitorPositions);
set(gcf,'color','w','name','The representation matrices')
numberOfSymmetry = size(symmetryIndex,2);
for i = 1:numberOfSymmetry
    subplot(2, 4, i);
    C = CC(i,:);
    n = sqrt(length(C));
    C = reshape(C,n,n);
    imagesc(C);
    axis equal;
    colormap(b2r(min(min(CC)),max(max(CC))))
    title(strcat('g_', num2str(i)))
end