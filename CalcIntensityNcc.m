%
%
function nccCoef = CalcIntensityNcc(d1Intensity1, d1Intensity2)
    nccCoef = 0;
    refLen = length(d1Intensity1); givenLen = length(d1Intensity2);
    if ( refLen == givenLen )
        nccCoef = ncc(d1Intensity1, d1Intensity2);
        return;
    end
	
    if ( refLen < givenLen )
        refMean = mean(d1Intensity1);
        nccCoef = ncc([d1Intensity1,repmat(refMean,1,givenLen-refLen)],d1Intensity2);
        return;
    end
    
    if ( refLen > givenLen )
        givenMean = mean(d1Intensity2);
        nccCoef = ncc(d1Intensity1,[d1Intensity2,repmat(givenMean,1,refLen-givenLen)]);
        return;
    end
end