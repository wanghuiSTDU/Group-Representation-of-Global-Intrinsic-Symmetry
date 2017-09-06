function [V,F] = read_off(filename)
%Read a triangular mesh from an OFF file.
%
%Input:
%filename: the filename of the triangular mesh in the OFF format.
%
%Output:
%V: an n*3 matrix specifying the position of the vertices, where n is the
%   number of vertices. 
%F: a m*3 matrix specifying the connectivity of the mesh, where m is the
%   number of triangular faces.
%
%Author: Hui Wang

%Open the file
fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

%Read the head of the file
str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end
str = fgets(fid);
sizes = sscanf(str, '%d %d %d', 3);
n = sizes(1);
m = sizes(2);

%Read the vertices
[A,cnt] = fscanf(fid, '%f %f %f', 3 * n);
if cnt ~= 3 * n
    warning('Problem in reading vertices.');
end
V = reshape(A, 3, cnt / 3);
V = V';

%Read the face indexes
[A, cnt] = fscanf(fid, '%d %d %d %d\n', 4 * m);
if cnt ~= 4 * m
    warning('Problem in reading faces.');
end
A = reshape(A, 4, cnt/4);
F = A(2:4,:)+1;
F = F'; 

%Close the file
fclose(fid);