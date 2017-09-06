function fun = indicatorFunction(V,F,index,para)
%Hui Wang

nIndex = length(index);
nVertex = size(V,1);
D = geodesicDistance_fastmarch_multiple(V,F,index);
for i = 1:nIndex
   D(i,i) = Inf;
end

dis = para * min(min(D));
fun = zeros(nVertex, nIndex);
for i = 1:nIndex
    d = geodesicDistance_fastmarch_single(V,F,index(i));
    dd = zeros(nVertex,1);
    dd(d < dis) = 1 \ sum(d < dis); %In practice, the indicator functions are discretized as density functions. 
    fun(:,i) = dd;
end