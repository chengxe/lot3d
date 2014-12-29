function [desireCoefs] = CalcEpipolarLineCoefs(stStereoModel, fmi, actualPt )
%
% actualPt - 2-by-1 Image point in actual image
%
    desireCoefs = stStereoModel.fundamentals{fmi} * [reshape(actualPt, 2, 1); 1];
    
end