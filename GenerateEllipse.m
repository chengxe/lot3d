function ellipse = GenerateEllipse(cov)
% generate an ellipse structure
%
%   input:
%       cov - covariance matrix of group of pixels
%
%   output:
%       ellipse - an ellipse structure, with components:
%               .radii - vector of radius of semi axis, ascend
%               .axis - 2-by-2, columns of semi axis, corr. to radii
%
    [ve, va] = eig(cov);
    
    [vaSorted, idxOrg] = sort(diag(va), 'ascend');
    veSorted = ve(:,idxOrg);
    
    ellipse.radii = sqrt(vaSorted);
    ellipse.axis = veSorted;
    ellipse.theta = atan2(veSorted(2,2),veSorted(1,2));
end