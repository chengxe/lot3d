function [measurements] = GetMeasurement(inImage)

    inImage = im2double(inImage); bw = inImage>0.5;

    [~, measurements] = InlineCalcRegionCentroid(bw, 20);
    
    pixelsWeightsOnIntensity = ones(size(inImage));
    for r = 1 : length(measurements)
        measurements(r).weights = pixelsWeightsOnIntensity(sub2ind(size(inImage),measurements(r).pixels(2,:),measurements(r).pixels(1,:)));
        measurements(r).mean = sum(measurements(r).pixels .* repmat(measurements(r).weights, 2, 1), 2) / sum(measurements(r).weights);
        measurements(r).cov  = ( ((measurements(r).pixels - repmat(measurements(r).mean,1,size(measurements(r).pixels,2))) .* repmat(measurements(r).weights, 2, 1)) * ...
                                 (measurements(r).pixels - repmat(measurements(r).mean,1,size(measurements(r).pixels,2)))' ) / sum(measurements(r).weights);
        measurements(r).cov = measurements(r).cov*4;
        measurements(r).ellipse = InlineGenerateEllipse(measurements(r).cov);   %可以用find(isnan(measurements(r).cov)) 或者 sum(sum(isnan(measurements(r).cov))) > 0 做debug条件
        measurements(r).ellipse.center = measurements(r).mean;
        measurements(r).ellipse.e = inv(measurements(r).cov);
        %measurements(r).pixel2ellipse = length(measurements(r).pixels) / (pi*measurements(r).ellipse.radii(1)*measurements(r).ellipse.radii(2));
        measurements(r).ellipse2pixel = (pi*measurements(r).ellipse.radii(1)*measurements(r).ellipse.radii(2)) / length(measurements(r).pixels);
    end
    
end


%-----------------------------------------
%
%-----------------------------------------
function [centroid, measurements] = InlineCalcRegionCentroid(bwImage, minArea)
    
    temp = bwlabel(bwImage);
    Sarea = regionprops(temp, 'Area');

    Sstats = regionprops(temp, 'PixelList');
    i = 1;
    for n=1:length(Sstats)
        if (Sarea(n).Area>minArea)
            measurements(i).area = Sarea(n).Area;
            measurements(i).pixels = Sstats(n).PixelList'; %第一列是列，第二列是行
            measurements(i).d1pixels = sub2ind(size(bwImage), measurements(i).pixels(2,:), measurements(i).pixels(1, :));
            measurements(i).d1limit = [max(measurements(i).d1pixels); min(measurements(i).d1pixels)];
            i = i+1;
        end
    end
    measurements = measurements';
    
   
    Sstats = regionprops(temp, 'centroid');
    centroid = [];
    i = 1;
    for n=1:length(Sstats)
        if (Sarea(n).Area>minArea)
            centroid = [centroid; Sstats(n).Centroid];
            measurements(i).centroid = Sstats(n).Centroid';
            i = i+1;
        end
    end
    
end


%
%
%
function ellipse = InlineGenerateEllipse(cov)
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