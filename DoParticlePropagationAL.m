function [tm0Particles, weights, particleMissed] = DoParticlePropagationAL(stStereoModel, stTrackingParameter, aTracker, tm0Frame, isDebug)
    
    tm0Particles = zeros(size(aTracker.particles)); weights = zeros(1, stTrackingParameter.numParticle); 
    
    if ( stTrackingParameter.noise.type == 1 )
        pNoise = stTrackingParameter.noise.sigma(1) * (0.5-rand(3,stTrackingParameter.numParticle));
    else
        pNoise = chol(stTrackingParameter.noise.sigma(2)) * randn(3,stTrackingParameter.numParticle);
    end
    
    prevTrajectory = aTracker.states(1:3, max(1,end-10):end);
    particleMissed = 0; 
    for n=1:stTrackingParameter.numParticle
        particles(n).associated = 0; particles(n).pal = -1; particles(n).pol = -1; 
        %--------------------------------------
        tm0Particles(:,n) =  stTrackingParameter.matTransition * aTracker.particles(:,n) + [pNoise(:,n); zeros(stTrackingParameter.lenState-3,1)];
        if ( IsOnWallOrOutCube(stTrackingParameter.cube, tm0Particles(1:3,n), stTrackingParameter.onoroutThreshold) )
            particleMissed = particleMissed + 1; weights(n) = 0; continue;
        end
        
        %-------------
        aSphere = bsxfun(@plus, stTrackingParameter.sphere, tm0Particles(1:3,n));
        okFlag = 1;
        for v = 1 : stTrackingParameter.numCamera
            d2Pixels = stStereoModel.cams(v).projection * [aSphere; ones(1,size(aSphere,2))];
            d2Pixels(1:2, :) = d2Pixels(1:2, :) ./ repmat(d2Pixels(3, :), 2, 1); d2Pixels = floor(d2Pixels(1:2,:)); d2Pixels = unique(d2Pixels','rows')';
            
            [occupiedRgns{v}, overlappedRatios{v}] = InlineFindOccupiedRegion(size(tm0Frame.cams(v).image),tm0Frame.cams(v).regions, d2Pixels);
            if (isempty(occupiedRgns{v}))
                particleMissed = particleMissed + 1; weights(n) = 0; okFlag = 0; break;
            end     
        end
        if ( 0 == okFlag ) continue; end

        for v = 1 : stTrackingParameter.numCamera
            particles(n).associated = 1; particles(n).camsMeasurements{v} = [];
            
            for idx = 1 : length(occupiedRgns{v})
                region = tm0Frame.cams(v).regions(occupiedRgns{v}(idx));
                
                particles(n).camsMeasurements{v} = [particles(n).camsMeasurements{v}, occupiedRgns{v}(idx)];
                %--------------------------------------------------------------------
                pal = overlappedRatios{v}(idx) * ...
                      InlineCalcIntensityNcc(aTracker.referee.cams(v).intensities, GetPixelIntensity(tm0Frame.cams(v).image, region.pixels));
                  
                if ( stTrackingParameter.haso )
                    pol = AttitudeLikelihood(stStereoModel.cams(v), tm0Frame.cams(v).image, region.ellipse, aTracker.referee.cams(v).ellipse, 0);
                else
                    pol = 20;
                end

                particles(n).weight.cams{v}(idx,:) = [pal, pol];
            end
        end
        
    end
    %------------------------------------------------------------------------------
    
    idxesAssociated = find(cat(2,particles.associated)==1); 
    if ( isempty(idxesAssociated) )
        particleMissed = stTrackingParameter.numParticle*2; weights = 1e-5; return; end
    for v = 1 : stTrackingParameter.numCamera
        %rgnsAssociated{v} = cat(2, particles(idxesAssociated).camsMeasurements{v}); rgnsAssociated{v} = unique(rgnsAssociated{v});
        rgnsAssociated{v} = [];
        for idx = idxesAssociated
            rgnsAssociated{v} = [rgnsAssociated{v}, particles(idx).camsMeasurements{v}]; 
        end
        rgnsAssociated{v} = unique(rgnsAssociated{v});
    end
    
    %-----------------------------------------
    for i=1:length(idxesAssociated)
        aParticle = particles(idxesAssociated(i));
        for v = 1 : stTrackingParameter.numCamera
            aParticle.camsMeasurements{v} = unique(aParticle.camsMeasurements{v}); %一个actual对应多个desire,会在actual产生重复
            for r=1:length(rgnsAssociated{v})
                idx = find(aParticle.camsMeasurements{v}==rgnsAssociated{v}(r));
                if ( ~isempty(idx) )
                    camsWeight{v}.pal(i,r) = aParticle.weight.cams{v}(idx,1);
                    camsWeight{v}.pol(i,r) = aParticle.weight.cams{v}(idx,2);
                else
                    camsWeight{v}.pal(i,r) = 0;
                    camsWeight{v}.pol(i,r) = 100;
                end
            end
        end
    end
    
    %-----------------------------------------
    for v = 1 : stTrackingParameter.numCamera
        
        [~, idxCamsRgn{v}] = sort(sum(camsWeight{v}.pal), 'descend');
    end
    
    %-----------------------------------------
    for i=1:length(idxesAssociated)
        aParticle = particles(idxesAssociated(i));
        
        for v = 1 : stTrackingParameter.numCamera
            idx = find(aParticle.camsMeasurements{v}==rgnsAssociated{v}(idxCamsRgn{v}(1)));
            if ( ~isempty(idx) )
                particles(idxesAssociated(i)).pal = 2*aParticle.weight.cams{v}(idx,1);
                particles(idxesAssociated(i)).pol = 0.5*aParticle.weight.cams{v}(idx,2);
            else
                idx = ismember(aParticle.camsMeasurements{v}, rgnsAssociated{v}(idxCamsRgn{v}(2:end)));
                idx = find(idx); 
                particles(idxesAssociated(i)).pal = aParticle.weight.cams{v}(idx(1),1);
                particles(idxesAssociated(i)).pol = aParticle.weight.cams{v}(idx(1),2);
            end            
        end
    end
    
    %-----------------------------------------    
    pals = cat(2, particles.pal); pols = cat(2, particles.pol);
    w1 = pals(idxesAssociated);
    if max(pols(idxesAssociated)) ~= min(pols(idxesAssociated))
        w2 = -(pols(idxesAssociated)-max(pols(idxesAssociated))) / (max(pols(idxesAssociated))-min(pols(idxesAssociated)));
    else
        w2 = pols(idxesAssociated) / max(pols(idxesAssociated));
    end
    
    %--------------------------------------------
    weights(idxesAssociated) = w1 .* exp(3*w2);

end

%
%
%
function [occupiedIdxes, overlappedRatio] = InlineFindOccupiedRegion(imgSize, regions, occupiedPoints)

    occupiedIdxes = []; overlappedRatio = [];
    if ( ~isempty(find(occupiedPoints(1,:)>imgSize(1))) || ~isempty(find(occupiedPoints(2,:)>imgSize(2))) || ...
         ~isempty(find(occupiedPoints(1,:)<1)) || ~isempty(find(occupiedPoints(2,:)<1)) )
        return;
    end
    try
        occupiedPixelList = sub2ind(imgSize, occupiedPoints(2,:), occupiedPoints(1, :));
        regionLimit = cat(2, regions.d1limit);
        idx1 = find(regionLimit(1,:)>=min(occupiedPixelList)); idx2 = find(regionLimit(2,:)<=max(occupiedPixelList));
        idxes = intersect(idx1, idx2);
        for i=1:length(idxes)
            overlapped = regions(idxes(i)).d1pixels(ismember(regions(idxes(i)).d1pixels, occupiedPixelList));
            if ( 0 == isempty(overlapped) )
                occupiedIdxes = [occupiedIdxes, idxes(i)];
                overlappedRatio(length(occupiedIdxes)) = length(overlapped) / length(regions(idxes(i)).d1pixels);
            end
            clear overlapped;
        end
    catch
        warning('max [%d;%d], min [%d;%d]',max(occupiedPoints(1,:)), max(occupiedPoints(2,:)), min(occupiedPoints(1,:)), min(occupiedPoints(2,:)));
        occupiedIdxes = []; overlappedRatio = []; return;
    end
end

%
%
%
function nccCoef = InlineCalcIntensityNcc(refIntensities, givenIntensities)
    nccCoef = 0;
    refLen = length(refIntensities); givenLen = length(givenIntensities);
    if ( refLen == givenLen )
        nccCoef = ncc(refIntensities, givenIntensities);
        return;
    end
	
    if ( refLen < givenLen )
        refMean = mean(refIntensities);
        nccCoef = ncc([refIntensities,repmat(refMean,1,givenLen-refLen)],givenIntensities);
        return;
    end
    
    if ( refLen > givenLen )
        givenMean = mean(givenIntensities);
        nccCoef = ncc(refIntensities,[givenIntensities,repmat(givenMean,1,refLen-givenLen)]);
        return;
    end
end
