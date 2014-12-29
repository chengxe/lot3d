function stBackgroundModel = GBMInitial(imageSeries, totalImage)
% stBackgroundModel.size = [row, col, window, total]
%
    stBackgroundModel = struct('size', [], ...
                               'mask', [], ...
                               'window', [], ...
                               'oIdx', [], ...
                               'sum', [], ...
                               'mean', [], ...
                               'std', []);
                               
    fprintf(1, 'Initialize background model...');
    
    stBackgroundModel.size = [size(imageSeries), totalImage];
    %stBackgroundModel.mask = zeros(stBackgroundModel.size(1:2));
    %stBackgroundModel.mask(50:stBackgroundModel.size(1)-50, 50:stBackgroundModel.size(2)-50) = 1;
    stBackgroundModel.mask = ones(stBackgroundModel.size(1:2));
    
    stBackgroundModel.window = imageSeries .* repmat(stBackgroundModel.mask, [1, 1, size(imageSeries,3)]);
    
    stBackgroundModel.sum = zeros(stBackgroundModel.size(1:2));
    for i=1:stBackgroundModel.size(3)
        stBackgroundModel.sum  = stBackgroundModel.sum + stBackgroundModel.window(:,:, i);
    end
    stBackgroundModel.oIdx = 1;
    
    stBackgroundModel.mean = stBackgroundModel.sum / stBackgroundModel.size(3);
    stBackgroundModel.std = zeros(stBackgroundModel.size(1:2));
    for i=1:stBackgroundModel.size(3);
        stBackgroundModel.std = stBackgroundModel.std + (stBackgroundModel.window(:,:,i) - stBackgroundModel.mean) .^ 2;
    end
    stBackgroundModel.std = (stBackgroundModel.std / (stBackgroundModel.size(3)-1)) .^ 0.5;
    stBackgroundModel.std(find(stBackgroundModel.std == 0)) = 1;   %%%%%%%%%%%%
    
    
    fprintf(1, 'done.\n');
end