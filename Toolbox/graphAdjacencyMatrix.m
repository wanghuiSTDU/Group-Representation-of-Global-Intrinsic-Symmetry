function W = graphAdjacencyMatrix(F)
%Compute the graph adjacency matrix of a triangular mesh.
%
%Input: 
%F: a m*3 matrix specifying the connectivity of the mesh, where m is the
%   number of triangular faces.
%
%Output:
%W: the n*n graph adjacency matrix. 
%
%Author: Hui Wang

n = max(max(F));
m = size(F,1);
W = sparse(n,n);

for ii = 1:m
   i = F(ii,1);   j = F(ii,2);   k = F(ii,3);
   
   W(i,j) = 1.0;
   W(i,k) = 1.0;

   W(j,i) = 1.0;
   W(j,k) = 1.0;

   W(k,j) = 1.0;
   W(k,i) = 1.0;
end