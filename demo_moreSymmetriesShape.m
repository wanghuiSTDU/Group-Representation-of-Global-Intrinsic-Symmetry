%Find all global intrinsic symmetries based on functional map for a table
% Hui Wang

clear;
clc;
close all;

%Add path
addpath(genpath(pwd));

%Read mesh
[V,F] = read_off('table.off');
fprintf('Loading triangular mesh with %d vertices and %d triangles.\n', size(V,1), size(F,1));

%Compute eigenvalues and eigenvectors of the Laplace-Beltrami operator 
fprintf('Eigen-decomposition of the Laplace-Beltrami operator.\n');
L = laplacian_matrix_cotan(V, F);
[eigvector, eigvalue] = eigs(L, 100, 'sm');
eigvalue = diag(eigvalue);
[eigvalue, index] = sort(eigvalue);
eigvector = eigvector(:, index);

%Repeated eigenvalues
fprintf('Computation of the voted matrices.\n');
reEigval{1} = 1;
reEigval{2} = [2,3];
reEigval{3} = [4,5];
reEigval{4} = [6,7];
reEigval{5} = [8];
reEigval{6} = [9,10];
reEigval{7} = [11,12];
reEigval{8} = [13,14];
reEigval{9} = [15,16];
reEigval{10} = [17];
reEigval{11} = [18,19];
reEigval{12} = [20];

%symmetry orbit
index = [61     63     4611     4645      9186     9187     13638    13649];

D = geodesicDistance_fastmarch_multiple(V,F,index);
[p1,p2] = generatePairsNew(8,2);
[pairs1,pairs2] = consistentPairs(p1,p2,D,0.9);

%Constraints functions
funIndicator = indicatorFunction(V,F,index,0.1);
hks = HKS(eigvector, eigvalue);
wks = compute_WKS(eigvector, eigvalue, 100);

%Functions to coefficents
numberOfEigenvector = 20;
eigvectorUsed = eigvector(:, 1:numberOfEigenvector);
funIndicatorC = eigvectorUsed' * funIndicator;
hksC = eigvectorUsed' * hks;
wksC = eigvectorUsed' * wks;

CC = [];
number = 1;
for ii = 1:size(pairs1,1)
   index1 = pairs1(ii,:);
   index2 = pairs2(ii,:);
   
   %Compute the matrix C
    C = zeros(numberOfEigenvector,numberOfEigenvector);
    paraHKS_WKS = 0.01;
    for i = 1:12
        in = reEigval{i};
        A1 = funIndicatorC(in, index1);
        B1 = funIndicatorC(in, index2);
        hksCi = hksC(in, :);
        wksCi = wksC(in, :);
        
        A = [A1, paraHKS_WKS * hksCi, paraHKS_WKS * wksCi];
        B = [B1, paraHKS_WKS * hksCi, paraHKS_WKS * wksCi];
        
        [u,s,v] = svd(B * A');
        R = u * v';
        C(in,in) = R;
    end
    
    CC(number,:) = C(:);
    number = number + 1;
end

%Mean-shift clustering
fprintf('Mean-shift clustering.\n');
h = 0.75 * mean(sqrt(sum(CC.^2,2)));
[clustCent,data2cluster,cluster2dataCell] = MeanShiftCluster(CC',h);
numberOfCluster = length(cluster2dataCell);
for i = 1:numberOfCluster
  numberOfEachCluster(i) = length(cluster2dataCell{i});
end
[v,in] = sort(numberOfEachCluster,'descend');

for i = 1:8
    clusterIndex = cluster2dataCell{in(i)};
    centerC = clustCent(:,in(i));
    D_C = zeros(length(clusterIndex),1);
    for j = 1:length(clusterIndex)
        C_j = CC(clusterIndex(j),:)';
        D_C(j) = norm(centerC - C_j, 'fro');
    end
    usedIndex = find(D_C == min(D_C));
    C_used(i,:) = CC(clusterIndex(usedIndex(1)),:);
end

%Refine by ICP
fprintf('ICP refinement.\n');
atria = nn_prepare(eigvectorUsed);
for j = 1:8
   C0 = reshape(C_used(j,:), numberOfEigenvector, numberOfEigenvector);
   [C, symIdx(:,j)] = refine_C(atria, C0, eigvectorUsed, 20, 10);    
   [C_usedRefined(j,:)] = reshape(C,1, numberOfEigenvector^2); 
end

%Visulization
fprintf('Visualization.\n');
load tableColor.txt
%Indices of the symmetries here don't coincide with the used ones in the paper
plot_moreSymmetriesShape(V, F, tableColor, index, C_usedRefined, symIdx);
fprintf('Done.\n');
