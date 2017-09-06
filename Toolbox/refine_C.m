function [C, symIdx] = refine_C(atria, C0, eigvector, nbr_iter, samplingInterval)
%The Iterative Closest Point (ICP) process.
%
%Input:
%atria: the data structure for nearest neighbors searching
%C0: the original symmetry representation matrix.
%eigvector: the eigenvectors of the Laplace-Beltrami operator.
%nbr_iter: the number of iterations, the default value is 20.
%samplingInterval: the sampling interval for ICP refinement, the default
%                  value is 10. 
%
%Output:
%C: the final symmetry representation matrix after ICP refinement. 
%symIdx: indexes of the point-to-point correspondence of the symmetry. 
%
%Author: Hui Wang

%Input parameters
if nargin < 3
    disp('Not enough parameters!');
elseif nargin < 4
    nbr_iter = 20;
    samplingInterval = 10;
elseif nargin < 5
    samplingInterval = 10;
elseif nargin >= 6
    disp('Too many parameters!');
end

%Initializations
queryPoints = eigvector(1:samplingInterval:end,:);
C = C0;

%ICP refinements
for i = 1:nbr_iter
    q = (C * queryPoints')';
    [symIdx0, dis] = nn_search(eigvector, atria, q, 1);
    imagePoints = eigvector(symIdx0,:);
    [U,S,V] = svd(imagePoints' * queryPoints);
    C = U * V';
end

%Final results
[symIdx, dis] = nn_search(eigvector, atria, (C * eigvector')', 1);  