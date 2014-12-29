function intensities = GetPixelIntensity(inImage, pixels)
%get the pixels intensities from image
%
%   input:
%       inImage - input original image
%       pixels - 2-by-n pixels list, first row is column(x axis), second is row(y axis).
%   
%   output:
%       intensities - n dims vector of intensities.
%
    pixels = floor(pixels);

    [Ry, Cx] = size(inImage);
    if ( (max(pixels(1,:)) > Cx || min(pixels(1,:)) < 1) || ...
         (max(pixels(2,:)) > Ry || min(pixels(2,:)) < 1) )   
        intensities = zeros([Ry, Cx]);
    else
        pixelList = sub2ind([Ry, Cx], pixels(2,:), pixels(1, :));
        inImage = double(inImage);
        intensities = inImage(pixelList);
    end
end