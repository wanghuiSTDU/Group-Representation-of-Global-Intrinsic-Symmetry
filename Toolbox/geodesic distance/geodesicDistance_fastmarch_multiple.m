function D = geodesicDistance_fastmarch_multiple(V,F,point_list)

%%Construct the geodesic distance from vertex point_list (a vector of index of vertex), the geodesic distance is
%approximated by the shortest distance on the mesh graph
%Input:  
%V ---the coordinates of the vertices of the triangle mesh
%F ---the fact of the triangle mesh
%point_list ---the index of the vertex
%Output: d ---the geodesic distance from the vertex i, d(i,j) is
%the geodesic distance between vertex point_list(i) and vertex
%point_list(j)
%%Example:
% [V F] = mRead();
% d = geodesicDistance_fastmarch_multiple(V,F,1:5);

%Hui Wang, wanghui19841109@gmail.com, First wirtten on April, 28, 2011


D = fastmarch(F,V(:,1),V(:,2),V(:,3),point_list,'multiple'); %This is copyed forn Alex Bronstein's code.