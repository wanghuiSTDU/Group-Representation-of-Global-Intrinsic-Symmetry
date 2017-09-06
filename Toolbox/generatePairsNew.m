%Hui Wang, Feb. 2, 2017
function [pairs1,pairs2] = generatePairsNew(n,m)
%Input:
%n: is the number of points
%m: is the number of pairs (n >= m)

pairs1 = combntns(1:n,m);
pairs2 = [];
for i = 1:size(pairs1,1)
   c = perms(pairs1(i,:));
   pairs2 = [pairs2;c];
end