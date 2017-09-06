function L = laplacian_matrix_cotan(V, F)
%Compute the Laplacian matrix using the classical cotangent weight scheme
%without area normalization.
%
%Input:
%V: an n*3 matrix specifying the position of the vertices, where n is the
%   number of vertices. 
%F: a m*3 matrix specifying the connectivity of the mesh, where m is the
%   number of triangular faces.
%
%Output:
%L: the n*n Laplacian matrix. 
%
%Author: Hui Wang

n = size(V,1);
m = size(F,1);
L = sparse(n,n);

for ii = 1:m
   i = F(ii,1);   j = F(ii,2);   k = F(ii,3);
   
   vi = V(i,:);   vj = V(j,:);   vk = V(k,:);
   
   alpha = myangle(vk - vi,vk - vj);
   beta = myangle(vj - vi,vj - vk);
   gama = myangle(vi - vj,vi - vk);

   L(i,j) = L(i,j) + cot(alpha);
   L(i,k) = L(i,k) + cot(beta);

   L(j,i) = L(j,i) + cot(alpha);
   L(j,k) = L(j,k) + cot(gama);

   L(k,j) = L(k,j) + cot(gama);
   L(k,i) = L(k,i) + cot(beta);
end    

L = 0.5 * (diag(sum(L,2)) - L);


function beta = myangle(u,v)

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );