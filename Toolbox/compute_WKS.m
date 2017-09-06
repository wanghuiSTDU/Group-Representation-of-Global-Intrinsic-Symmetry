% This is the implementation of the Wave Kernel Signature described
% in the paper:
% 
%    The Wave Kernel Signature: A Quantum Mechanical Approach To Shape Analysis 
%    M. Aubry, U. Schlickewei, D. Cremers
%    In IEEE International Conference on Computer Vision (ICCV) - Workshop on 
%    Dynamic Shape Capture and Analysis (4DMOD), 2011
% 
% Please refer to the publication above if you use this software. 


function [WKS] = compute_WKS(PHI, E, nWKS, wks_variance)
%% 
% compute the Wave Kernel Signature of triangle mesh given by  [vertices,faces]
%   
%   INPUT:
%   PHI is the (number of vertices x 300) matrix of LB eigenfunctions 
%   E is the vector of LB eigenvalues (by default of size 300 x 1)
%   nWKS is the number of Wave Kernel Singature (by default is 100)
%   wks_variance is variance of the WKS gaussian (wih respect to the difference of the two first eigenvalues)
%   
%   OUTPUT:
%   WKS is the (number of vertices) x 100 WKS matrix 
%
%
%   The main parameter to adjust depending on your task is wks_variance


%% parameters
if nargin < 3
    nWKS = 100; % number of evaluations of WKS
    wks_variance = 6; % variance of the WKS gaussian (wih respect to the 
    % difference of the two first eigenvalues). For easy or precision tasks 
    % (eg. matching with only isometric deformations) you can take it smaller
end
if nargin < 4
    wks_variance = 6; 
end

%% compute WKS 
num_vertices = size(PHI,1);
nWKS = min (size(PHI,2), nWKS);
WKS=zeros(num_vertices,nWKS);

log_E=log(max(abs(E),1e-6))';
e=linspace(log_E(2),(max(log_E))/1.02,nWKS);  
sigma=(e(2)-e(1))*wks_variance;

C = zeros(1,nWKS); %weights used for the normalization of f_E

for i = 1:nWKS
    WKS(:,i) = sum(PHI.^2.*...
        repmat( exp((-(e(i) - log_E).^2) ./ (2*sigma.^2)),num_vertices,1),2);
    C(i) = sum(exp((-(e(i)-log_E).^2)/(2*sigma.^2)));
end

% normalize WKS
WKS(:,:) = WKS(:,:)./repmat(C,num_vertices,1);
