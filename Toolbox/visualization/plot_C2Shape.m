function plot_C2Shape(V,F,eigvector,C,symPairs,sampleIndex, sampleIndexSym)
%Hui Wang

a=get(0);
figure('position',a.MonitorPositions);

subplot(1,4,1);
visualize_symmetry_lines(V,F,symPairs([1;3]),symPairs([2;4]));
view([-22 -6]);
title('Symmetry pairs');

subplot(1,4,2);
visualize_symmetry_lines(V,F,sampleIndex,sampleIndexSym);
title('Point to point correspondence');
view([-22 -6]);

subplot(1,4,3);
Coeff = eigvector' * V(:,1);
orginalColor = eigvector * Coeff;
trimesh(F,V(:,1),V(:,2),V(:,3), orginalColor, 'FaceColor', 'interp', 'EdgeColor', 'none');
colormap jet;
axis off;
axis equal;
view([-22 -6]);
title('Original function')

subplot(1,4,4);
trasferColor = eigvector * C * Coeff;
trimesh(F,V(:,1),V(:,2),V(:,3), trasferColor, 'FaceColor', 'interp', 'EdgeColor', 'none');
colormap jet;
axis off;
axis equal;
view([-22 -6]);
title('Symmetry transfer function')

figure;
imagesc(C)
set(gcf,'color','w')
axis equal;
colormap(b2r(min(min(C)),max(max(C))))
title('Representation matrix')

