function regions = ReCalcRegions(regions, inImage, stBackModel)
%calculate others components of regoins
%
%   input:
%       regions - n 'region' structure, each with components:
%                   .area
%                   .pixels
%                   .centroid
%       inImage - original gray image
%       stBackModel - background model, with component:
%                       .mask - mask of edge at each side
%                       .mean - mean intensity of each pixel
%                       .std - std of intensity at each pixel
%   output:
%       regions - n 'region' structure, each with components:
%                   .area
%                   .pixels
%                   .centroid
%                   .mean
%                   .cov
%                   .ellipse.center - center of an ellipse
%                   .ellipse.radii - vector of radius of semi axis, ascend
%                   .ellipse.axis - 2-by-2, columns of semi axis, corr. to radii
%                   .ellipse.e - inverse of covariance
%                   .pixel2ellipse - metric for how an ellipse fitting the region pixels
%

    inImage = im2double(inImage); inImage = inImage .* stBackModel.mask;
    
     %------------------------------------------------------------------------
     pixelWeights = abs((inImage-stBackModel.mean) ./ stBackModel.std);  emptyRgns = []; %weights = [];
    for r=1:length(regions)
        regions(r).weights = pixelWeights(sub2ind(size(inImage),regions(r).pixels(2,:),regions(r).pixels(1,:)));
        if ( max(regions(r).weights) < 3 )
            emptyRgns = [emptyRgns, r]; continue;
        end
%         weights = [weights, max(regions(r).weights)];
        regions(r).mean = sum(regions(r).pixels .* repmat(regions(r).weights, 2, 1), 2) / sum(regions(r).weights);
        regions(r).cov  = ( ((regions(r).pixels - repmat(regions(r).mean,1,size(regions(r).pixels,2))) .* repmat(regions(r).weights, 2, 1)) * ...
                                 (regions(r).pixels - repmat(regions(r).mean,1,size(regions(r).pixels,2)))' ) / sum(regions(r).weights);
        regions(r).cov = regions(r).cov * 5;
        regions(r).ellipse = GenerateEllipse(regions(r).cov);   %可以用find(isnan(regions(r).cov)) 或者 sum(sum(isnan(regions(r).cov))) > 0 做debug条件
        regions(r).ellipse.center = regions(r).mean;
        regions(r).ellipse.e = inv(regions(r).cov);
        %regions(r).pixel2ellipse = length(regions(r).pixels) / (pi*regions(r).ellipse.radii(1)*regions(r).ellipse.radii(2));
        regions(r).ellipse2pixel = (pi*regions(r).ellipse.radii(1)*regions(r).ellipse.radii(2)) / length(regions(r).pixels);
    end    
    regions(emptyRgns) = []; %weights = [];
%     %------------------------------------------------------------------------
%     p2e = cat(1,regions.pixel2ellipse);
%     reThresholdIdx = find(p2e<(mean(p2e)-2*std(p2e)));
%     for r=1:length(reThresholdIdx)
%         newRegions = ReThresholdRegion(inImage, regions(reThresholdIdx(r)).pixels);
%         if (isempty(newRegions)) 
%             continue; end
%         emptyRgns = [];
%         for rr=1:length(newRegions)
%             newRegions(rr).weights = pixelWeights(sub2ind(size(inImage),newRegions(rr).pixels(2,:),newRegions(rr).pixels(1,:)));
%             if ( max(newRegions(rr).weights) < 5 )
%                 emptyRgns = [emptyRgns, rr]; continue;
%             end
% %             weights = [weights, max(newRegions(rr).weights)];
%             newRegions(rr).mean = sum(newRegions(rr).pixels .* repmat(newRegions(rr).weights, 2, 1), 2) / sum(newRegions(rr).weights);
%             newRegions(rr).cov  = ( ((newRegions(rr).pixels - repmat(newRegions(rr).mean,1,size(newRegions(rr).pixels,2))) .* repmat(newRegions(rr).weights, 2, 1)) * ...
%                                      (newRegions(rr).pixels - repmat(newRegions(rr).mean,1,size(newRegions(rr).pixels,2)))' ) / sum(newRegions(rr).weights);
%             newRegions(rr).cov = newRegions(rr).cov * 5;
%             newRegions(rr).ellipse = GenerateEllipse(newRegions(rr).cov);
%             newRegions(rr).ellipse.center = newRegions(rr).mean;
%             newRegions(rr).ellipse.e = inv(newRegions(rr).cov);
%             newRegions(rr).pixel2ellipse = length(newRegions(rr).pixels) / (pi*newRegions(rr).ellipse.radii(1)*newRegions(rr).ellipse.radii(2));
%         end
%         newRegions(emptyRgns) = [];
%         if ( isempty(newRegions) )
%             continue; end
%         if ( length(newRegions) == 1 )
%             regions(reThresholdIdx(r)) = newRegions(1);
%         else
%             regions(reThresholdIdx(r)) = newRegions(1);
%             regions(length(regions)+1:length(regions)+1+length(newRegions)-1) = newRegions(2:length(newRegions));
%         end
%     end  

end

%
%
%
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

%
%
%
function [regions,outImage] = ReThresholdRegion(inImage, regionPixels)

    bwImage = zeros(size(inImage));
    
    ind1 = sub2ind(size(inImage), regionPixels(2,:), regionPixels(1,:));
    bwImage(ind1) = 1;
    
    se = strel('disk',5);
    bwImage = imdilate(bwImage, se);
    
    regions = CalcRegionCentroid(bwImage, 20);
    if ( isempty(regions) )
        outImage = inImage; end
    [regions, outImage] = RegionRefined(inImage, regions, 1.7);
%     for err=1.6:0.2:2
%         [outRegions, outImage] = RegionRefined(inImage, regions, err);
%         
%     end
    
end

%-----------------------------------------
%
%-----------------------------------------
function [regions, centroid] = CalcRegionCentroid(bwImage, minArea)
%
% inImage is an binary image.
% region. pixels
%         centroid
%
    
    temp = bwlabel(bwImage);
    Sarea = regionprops(temp, 'Area');

    Sstats = regionprops(temp, 'PixelList');
    i = 1;
    for n=1:length(Sstats)
        if (Sarea(n).Area>minArea)
            regions(i).area = Sarea(n).Area;
            regions(i).pixels = Sstats(n).PixelList'; %第一列是列，第二列是行
            i = i+1;
        end
    end
    if ( exist('regions','var') == 0 )
        regions = []; centroid = []; end
    regions = regions';
    
   
    Sstats = regionprops(temp, 'centroid');
    centroid = [];
    i = 1;
    for n=1:length(Sstats)
        if (Sarea(n).Area>minArea)
            centroid = [centroid; Sstats(n).Centroid];
            regions(i).centroid = Sstats(n).Centroid';
            i = i+1;
        end
    end
    
end

%-----------------------------------------
% 针对邻域内的原图像素灰度进行操作，
%-----------------------------------------
function [outRegions, outImage] = RegionRefined(inImage, inRegions, stdError)
% region(i).centroid
%          .area
%          .pixels
%
    
    outImage = zeros(size(inImage));
	for i=1:length(inRegions)

        ind1 = sub2ind(size(inImage), inRegions(i).pixels(2,:), inRegions(i).pixels(1,:));
        regionIntensities = inImage(ind1);
        regionMean = mean(regionIntensities);
        regionStd  = std(regionIntensities);
        
        ind2 = [];
        for j=1:length(ind1)
            if ( inImage(ind1(j)) <= regionMean-regionStd*stdError)
                ind2 = [ind2; ind1(j)];
            end
        end
        if ( length(ind2) > 5 )
            outImage(ind2) = 1;
        end
    end
    outRegions = CalcRegionCentroid(outImage, 5);
end