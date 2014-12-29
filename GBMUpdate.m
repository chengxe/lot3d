function outModel = GBMUpdate(inImage, inModel)
%
%
    outModel = inModel;
    
    outModel.sum = outModel.sum - outModel.window(:,:,outModel.oIdx);
    outModel.window(:,:,outModel.oIdx) = inImage .* outModel.mask;
    outModel.sum = outModel.sum + outModel.window(:,:,outModel.oIdx);
    %-------------------------------
    outModel.oIdx = outModel.oIdx +1;
    if ( outModel.oIdx > outModel.size(3) )
        outModel.oIdx = 1; end
    
    outModel.mean = outModel.sum / outModel.size(3);
    outModel.std = zeros(outModel.size(1:2));
    for i=1:outModel.size(3);
        outModel.std = outModel.std + abs(outModel.window(:,:,i) - outModel.mean);
    end
    outModel.std = 1.4826*(outModel.std / outModel.size(3));
    outModel.std(find(outModel.std == 0)) = 1;   %%%%%%%%%%%%
    
end