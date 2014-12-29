function [regions, bwImage] = GetObjectRegion(inImage, stBackModel)
% find connected regions in image, after background subtracion.
%
%   input:
%       inImage - an image
%       stBackModel - background model, with component:
%                       .mask - mask of edge at each side
%                       .mean - mean intensity of each pixel
%                       .std - std of intensity at each pixel
%
%   output:
%       regions - n 'region' structure, each with components:
%                   .area
%                   .pixels
%                   .centroid
%       bwImage - a binary image, by substracting background from input
%

     inImage = im2double(inImage); inImage = inImage .* stBackModel.mask;
     bwImage = abs((inImage-stBackModel.mean) ./ stBackModel.std) > 3; %用全局高斯背景是5
     
     se = strel('disk',5);
%      bwImage = imopen(bwImage, se);
     bwImage = imdilate(bwImage, se);
     
     [~, regions] = CalcRegionCentroid(bwImage, 10); %用全局高斯背景是15

     [regions, bwImage] = RegionRefined(inImage, regions, 1.7);
end

%-----------------------------------------
%
%-----------------------------------------
function [centroid, regions] = CalcRegionCentroid(bwImage, minArea)
% calculate region's centroid, area, and pixels
%
% input:
%       bwImage - input image is a binary image.
%       minArea - min area for a valid region
%
% output:
%       centroid - n-by-2 matrix of centroids of all the valid region
%       regions - n 'region' structure, each with components:
%                   .area
%                   .pixels
%                   .centroid
%

    
    temp = bwlabel(bwImage);
    Sarea = regionprops(temp, 'Area');

    Sstats = regionprops(temp, 'PixelList');
    i = 1;
    for n=1:length(Sstats)
        if (Sarea(n).Area>minArea)
            regions(i).area = Sarea(n).Area;
            regions(i).pixels = Sstats(n).PixelList'; %第一列是列，第二列是行
            regions(i).d1pixels = sub2ind(size(bwImage), regions(i).pixels(2,:), regions(i).pixels(1, :));
            regions(i).d1limit = [max(regions(i).d1pixels); min(regions(i).d1pixels)];
            i = i+1;
        end
    end
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
        if ( length(ind2) > 10 && length(ind2) < 200 )   % && length(ind2) < 100 这个条件是因为挥动布的影响
            outImage(ind2) = 1;
        end
    end
    [~, outRegions] = CalcRegionCentroid(outImage, 10);
end

