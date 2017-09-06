%Find the global intrinsic symmetries of a human model 
%Hui Wang
% 
clear;
clc;
close all;

%Add path
addpath('./Data');
addpath('./Toolbox');
addpath('./Toolbox/geodesic distance');
addpath('./Toolbox/nnsearch');
addpath('./Toolbox/visualization');

%Read mesh
[V, F] = read_off('horse.off');
fprintf('Loading triangular mesh with %d vertices and %d triangles.\n', size(V,1), size(F,1));

%Compute eigenvalues and eigenvectors of the Laplace-Beltrami operator 
fprintf('Eigen-decomposition of the Laplace-Beltrami operator.\n');
L = laplacian_matrix_cotan(V, F);
[eigvector, eigvalue] = eigs(L, 100, 'sm');
eigvalue = diag(eigvalue);
[eigvalue, index] = sort(eigvalue);
eigvector = eigvector(:, index);

%Repeated eigenvalues
reEigval{1} = [1];
reEigval{2} = [2];
reEigval{3} = [3 4 5];
reEigval{4} = [6];
reEigval{5} = [7];
reEigval{6} = [8];
reEigval{7} = [9];
reEigval{8} = [10];
reEigval{9} = [11 12];
reEigval{10} = [13];
reEigval{11} = [14 15];
reEigval{12} = [16];
reEigval{13} = [17];
reEigval{14} = [18];
reEigval{15} = [19 20];

%Heat Kernel Signatures, Wave Kernel Signatures, and pairs generation
fprintf('Computation of the constraint functions.\n');
[hks] = HKS(eigvector, eigvalue);
wks = compute_WKS(eigvector, eigvalue);
[maximaIndex,minimaIndex] = extrema(F, hks(:,end), 2);
extremaIndex = [maximaIndex;minimaIndex];
    
D = geodesicDistance_fastmarch_multiple(V,F,extremaIndex);
maxD = max(max(D));
    
DIS = dist(hks(extremaIndex,:)');
for i = 1:length(extremaIndex)
    DIS(i,i) = Inf;
end
    
[in1,in2,v] = find(DIS == min(min(DIS)));
 while 1
    if D(in1(1),in2(1)) > 0.3 * maxD
       DIS([in1,in2],:) = Inf;
       DIS(:,[in1,in2]) = Inf;
       break;
    else
       DIS(in1(1),in2(1)) = Inf;
       DIS(in2(1),in1(1)) = Inf;
       [in1,in2,v] = find(DIS == min(min(DIS)));
    end
 end
 d1 = D(in1(1),:);
 ii1 = find(d1 < 0.3 * maxD);
 DIS(ii1,:) = Inf;
 DIS(:,ii1) = Inf;
   
 d2 = D(in2(1),:);
 ii2 = find(d2 < 0.3 * maxD);
 DIS(ii2,:) = Inf;
 DIS(:,ii2) = Inf;
    
 [in3,in4,v] = find(DIS == min(min(DIS)));
 while 1
    if D(in3(1),in4(1)) > 0.3 * maxD
       DIS([in3,in4],:) = Inf;
       DIS(:,[in3,in4]) = Inf;
       break;
    else
       DIS(in3(1),in4(1)) = Inf;
       DIS(in3(1),in4(1)) = Inf;
       [in3,in4,v] = find(DIS == min(min(DIS)));
    end
end
   
inn = extremaIndex([in1;in3]);
DD = D([in1;in3], [in1;in3]);
for i = 1:4
    DD(i,i) = Inf;
end
minDD = min(min(DD));
para = 0.1;
%In practice, the indicator functions are discretized as density functions. 
d1 = geodesicDistance_fastmarch_single(V,F,inn(1));
f1 = 0 * d1;
f1(d1 < para * minDD) = 1 / sum(d1 < para * minDD);
    
d2 = geodesicDistance_fastmarch_single(V,F,inn(2));
f2 = 0 * d2;
f2(d2 < para * minDD) = 1 / sum(d2 < para * minDD);
    
d3 = geodesicDistance_fastmarch_single(V,F,inn(3));
f3 = 0 * d3;
f3(d3 < para * minDD) = 1 / sum(d3 < para * minDD);
    
d4 = geodesicDistance_fastmarch_single(V,F,inn(4));
f4 = 0 * d4;
f4(d4 < para * minDD) = 1 / sum(d4 < para * minDD);

numberOfEigenvector = 20;
eigvector = eigvector(:,1:numberOfEigenvector);
eigvalue = eigvalue(1:numberOfEigenvector);
hksC = eigvector' * hks;
wksC = eigvector' * wks;

%Band-by-band computation of the symmetry representation matrix
fprintf('Band-by-band computation of the representation matrix.\n');
C0 = zeros(numberOfEigenvector,numberOfEigenvector);
for i = 1:15
    in = reEigval{i};
    A1 = eigvector(:,in)' * [f1,f2,f3,f4];
    B1 = eigvector(:,in)' * [f2,f1,f4,f3];
        
    hksC_i = hksC(in,:);
    wksC_i = wksC(in,:);
    A = [A1, 0.01 * wksC_i, 0.01 * hksC_i];
    B = [B1, 0.01 * wksC_i, 0.01 * hksC_i];
        
    [u,s,v] = svd(B * A');
    R = u * v';
    C0(in,in) = R;
end

%Iterative Closest Point refinement
fprintf('ICP refinement.\n');
atria = nn_prepare(eigvector); 
[C, symIdx] = refine_C(atria, C0, eigvector, 20, 10);

%Visualization
fprintf('Visualization.\n');
plot_C2Shape(V,F,eigvector, C, inn',1:200:size(V,1), symIdx(1:200:size(V,1)));
fprintf('Done.\n');
