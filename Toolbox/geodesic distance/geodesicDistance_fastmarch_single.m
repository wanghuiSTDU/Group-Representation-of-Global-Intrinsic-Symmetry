function D = geodesicDistance_fastmarch_single(V,F,i)

%%Construct the geodesic distance from vertex point_list (a vector of index of vertex), the geodesic distance is
%approximated by the shortest distance on the mesh graph
%Input:  
%V ---the coordinates of the vertices of the triangle mesh
%F ---the fact of the triangle mesh
%i ---the index of the signal vertex
%Output: d ---the geodesic distance from the vertex i, d(i,j) is
%the geodesic distance between vertex i and vertex j
%%Example:
%[V F] = mRead();
%d = geodesicDistance_fastmarch_single(V,F,1);

%Hui Wang, wanghui19841109@gmail.com, First wirtten on April, 28, 2011

D0 = repmat(Inf,[size(V,1) 1]);
D0(i) = 0;
D = fastmarch(F,V(:,1),V(:,2),V(:,3),D0, struct('mode', 'single')); %This is copyed forn Alex Bronstein's code.