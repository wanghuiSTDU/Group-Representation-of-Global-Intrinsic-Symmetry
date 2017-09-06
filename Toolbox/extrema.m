function [maximaIndex,minimaIndex] = extrema(F, fun, k)
%Seeking the local maxima and minima points of the function, the maxima is
%defined as: whose value is larger than the values at its k-ring neighbor
%points.
%Input
%
%F: a m*3 matrix specifying the connectivity of the mesh, where m is the
%   number of triangular faces.
%fun: the functional values defined on the vertices.
%k: the number of the rings.
%
%Output
%maximaIndex: the indexes of the local maxima points.
%minmaIndex: the indexes of the local minima points.

%Author: Hui Wang

W = graphAdjacencyMatrix(F);
n = size(W,1);
if k > 1
  neighbourMatrix = spones(W^k) - speye(n);
elseif k == 1
  neighbourMatrix = W;
else
  error('k must be a negative integer!');
end

[row,col] = find(neighbourMatrix);
valPlus = fun(row) > fun(col);
flagMatrixPlus = sparse(row,col,valPlus);
valMinus = fun(row) < fun(col);
flagMatrixMinus = sparse(row,col,valMinus);

num = sum(neighbourMatrix,2);
numPlus = sum(flagMatrixPlus,2);
numMinus = sum(flagMatrixMinus,2);

maximaIndex = find(numPlus >= num);
minimaIndex = find(numMinus >= num);
