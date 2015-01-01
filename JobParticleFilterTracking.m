function updatedTrackers = JobParticleFilterTracking(stStereoModel, stTrackingParameter, trackers, idxesTracking, tm0Frame, t,isDebug, issimu)
%
%
%
%     rng('default');
    cntTrackers = length(trackers); updatedTrackers = [];
    for n=1:cntTrackers
%         disp(['tracker=',num2str(idxesTracking(n))]);
        %---------------------------------------------
        for retry = 1 : stTrackingParameter.threshold.retry
            aTracker = trackers(n);
            %---------------------------------------------
            aTracker.tm1 = aTracker.tm0; aTracker.tm0 = aTracker.tm0+1;
        
            [aTracker.particles, aTracker.weights, particleMissed] = DoParticlePropagationAL(stStereoModel, stTrackingParameter, aTracker, tm0Frame, 0);
            %----------------------------------------------------------------
            if (particleMissed > stTrackingParameter.threshold.missing.particle)
                aTracker.missing = aTracker.missing + 1;
            end
            if ( max(max(aTracker.weights)) <= 1e-4 )
                aTracker.missing = aTracker.missing + 1;
            end
            %-----------------------------
            if (aTracker.missing > stTrackingParameter.threshold.missing.target)
                aTracker.end = t-1; aTracker.activated = 0; 
                warning('tracker %d, time %d, retry %d',idxesTracking(n), t, retry);
                rng('shuffle'); continue;
            end
            aTracker.missing = 0;
            aTracker.weights = aTracker.weights ./ sum(aTracker.weights);
            newState = aTracker.particles * aTracker.weights'; newState(6:8) = aTracker.states(1:3,end); 
            %-------------------------------------------------------
            %是否在已经处理过的目标的重合距离之内
            %
            %-------------------------------------------------------
            if ( IsOnWallOrOutCube(stTrackingParameter.cube, newState(1:3), stTrackingParameter.onoroutThreshold) )
            %在外面
                aTracker.end = t-1; aTracker.activated = 0;  
                warning('tracker %d, time %d, retry %d',idxesTracking(n), t, retry);
                rng('shuffle'); continue;
            end
            
            %-------------------------------------------------------
            occupiedRgn = InlineFindExpectationOccupiedRgn(stStereoModel, stTrackingParameter, tm0Frame, newState);
            if ( isempty(occupiedRgn) )
            % 没有关联到量测 
                aTracker.end = t-1; aTracker.activated = 0;  
                warning('tracker %d, time %d, retry %d',idxesTracking(n), t, retry);
                rng('shuffle'); continue;
            end
            if ( stTrackingParameter.haso )
                %
                % 计算 orientation 和 candidate location
                %
                for v = 1 : stTrackingParameter.numCamera
                    gammas(v) = tm0Frame.cams(v).regions(occupiedRgn(v)).ellipse.radii(2) / tm0Frame.cams(v).regions(occupiedRgn(v)).ellipse.radii(1);
                    if ( tm0Frame.cams(v).regions(occupiedRgn(v)).ellipse2pixel > stTrackingParameter.threshold.merged )
                        gammas(v) = 1 / tm0Frame.cams(v).regions(occupiedRgn(v)).ellipse2pixel;
                    end
                end
                [~, ix] = sort(gammas, 'descend'); ix = ix(1:2); ix = sort(ix, 'ascend');
                d3Orientation = ReconstructOrientation(stStereoModel.cams(ix(1)), stStereoModel.cams(ix(2)),...
                    tm0Frame.cams(ix(1)).regions(occupiedRgn(ix(1))).ellipse,...
                    tm0Frame.cams(ix(2)).regions(occupiedRgn(ix(2))).ellipse, issimu);
                d3Location = BinocularReconstruction(stStereoModel.cams(ix(1)), stStereoModel.cams(ix(2)),...
                    tm0Frame.cams(ix(1)).regions(occupiedRgn(ix(1))).centroid,...
                    tm0Frame.cams(ix(2)).regions(occupiedRgn(ix(2))).centroid);

                if ( ~issimu )
                    [theta0, phi0, ~] = cart2sph(d3Orientation(1), d3Orientation(3), -d3Orientation(2));
                else
                    [theta0, phi0, ~] = cart2sph(d3Orientation(1), d3Orientation(2), d3Orientation(3));
                end
                newState(4) = theta0; newState(5) = phi0;
                %--------------------------
                dis = distance(d3Location, newState(1:3));
                if ( dis > stTrackingParameter.lenBody )
                   newState(1:3) = (d3Location+newState(1:3)) / 2;
                end
                %--------------------------
            end
            aTracker.states = [aTracker.states, newState];
            %----------------------------------------------------------------
            %[aTracker.particles, aTracker.weights] = InlineParticleResample(aTracker.particles, aTracker.weights, 'Systematic');
            aTracker.particles = repmat(aTracker.states(:,end),1,stTrackingParameter.numParticle);
            if ( stTrackingParameter.noise.type == 1 )
                aTracker.particles(1:3,:) = aTracker.particles(1:3,:) + stTrackingParameter.noise.sigma(1)*(0.5-rand(3,stTrackingParameter.numParticle));
            else
                aTracker.particles(1:3,:) = aTracker.particles(1:3,:) + chol(stTrackingParameter.noise.sigma(2))*randn(3,stTrackingParameter.numParticle);
            end            
            %----------------------------------------------------------------
            % 更新模板
            for v = 1 : stTrackingParameter.numCamera
                aTracker.referee.cams(v).ellipse = tm0Frame.cams(v).regions(occupiedRgn(v)).ellipse;
                aTracker.referee.cams(v).intensities = GetPixelIntensity(tm0Frame.cams(v).image, tm0Frame.cams(v).regions(occupiedRgn(v)).pixels);
            end

            aTracker.end = t;
            if ( aTracker.activated == 1 )
                break; end
            warning('tracker %d, time %d, retry %d',idxesTracking(n), t, retry);
            rng('shuffle');
        end
        updatedTrackers = [updatedTrackers; aTracker];
        clear aTracker; %toc
    end
    
end

%
%
%
function [newParticles, newWeights] = InlineParticleResample(inParticles, inWeights, strategy)

    numParticle = length(inWeights);
    switch strategy
        case 'Multinomial'
            idx = randsample(1:numParticle, numParticle, true, inWeights);
            
        case 'Systematic'
            cdf = min([0 cumsum(inWeights)], 1);
            cdf(end) = 1;
            
            u1 = rand / numParticle;
            [~, idx] = histc(u1 : 1/numParticle : 1, cdf);
    end

    newParticles = inParticles(:, idx);
    newWeights   = repmat(1/numParticle, 1, numParticle);
end

%
%
%
function [occupiedRgn] = InlineFindExpectationOccupiedRgn(stStereoModel, stTrackingParameter, tm0Frame, expectation)
    
    if ( stTrackingParameter.haso )
        [d3Orientation(1), d3Orientation(2), d3Orientation(3)] = sph2cart(expectation(4), expectation(5), 1);
        shape = GenerateBodyShape(stTrackingParameter.shape.profile, reshape(d3Orientation,3,1), expectation(1:3), 100);
        d3Upsilon = shape.points;
    else
        d3Upsilon = bsxfun(@plus, stTrackingParameter.sphere, expectation(1:3));
    end
    %----------------------------------------------------------------------------
    for v = 1 : stTrackingParameter.numCamera
        d2Pixels = stStereoModel.cams(v).projection * [d3Upsilon; ones(1,size(d3Upsilon, 2))];
        d2Pixels(1:2, :) = d2Pixels(1:2, :) ./ repmat(d2Pixels(3, :), 2, 1); d2Pixels = floor(d2Pixels(1:2,:)); d2Pixels = unique(d2Pixels','rows')';
        
        [occupiedRgns, overlappedSizes] = InlineFindOccupiedRegion(size(tm0Frame.cams(v).image),tm0Frame.cams(v).regions, d2Pixels);
        if (isempty(occupiedRgns))
            occupiedRgn = []; return;
        end
        [~, idx] = max(overlappedSizes);
        occupiedRgn(v) = occupiedRgns(idx);
    end

end

%
%
%
function [occupiedIdxes, overlappedSizes] = InlineFindOccupiedRegion(imgSize, regions, occupiedPoints)
    
    occupiedIdxes = []; overlappedSizes = [];
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
                overlappedSizes = [overlappedSizes, length(overlapped)];
            end
            clear overlapped;
        end
    catch
        warning('max [%d;%d], min [%d;%d]',max(occupiedPoints(1,:)), max(occupiedPoints(2,:)), min(occupiedPoints(1,:)), min(occupiedPoints(2,:)));
        occupiedIdxes = []; overlappedSizes = []; return;
    end

end

%
%
%
function DisplaySubImage(stStereoModel, actualImage, desireImage, curState, type)
    scale = 10; subSize = 25;
    
    curActualPixel = stStereoModel.matActualProjection * [curState(1:3);1]; curActualPixel = floor(curActualPixel /curActualPixel(3));
    curDesirePixel = stStereoModel.matDesireProjection * [curState(1:3);1]; curDesirePixel = floor(curDesirePixel /curDesirePixel(3));
    
    [Ry, Cx] = size(actualImage);
    if ( (max(curActualPixel(1,:)) > Cx-subSize || min(curActualPixel(1,:)) < subSize+1) || ...
         (max(curActualPixel(2,:)) > Ry-subSize || min(curActualPixel(2,:)) < subSize+1) )   
        return;
    end
    if ( (max(curDesirePixel(1,:)) > Cx-subSize || min(curDesirePixel(1,:)) < subSize+1) || ...
         (max(curDesirePixel(2,:)) > Ry-subSize || min(curDesirePixel(2,:)) < subSize+1) )   
        return;
    end
     
    try
        subActualImage = actualImage(curActualPixel(2)-subSize:curActualPixel(2)+subSize,curActualPixel(1)-subSize:curActualPixel(1)+subSize,:);
        subDesireImage = desireImage(curDesirePixel(2)-subSize:curDesirePixel(2)+subSize,curDesirePixel(1)-subSize:curDesirePixel(1)+subSize,:);
    catch
        warning('display sub-image error: state=(%s)',num2str(reshape(curState,1,10)));
        return;
    end
    subActualImage = imresize(subActualImage, scale);
    subDesireImage = imresize(subDesireImage, scale);
    
    [x,y,z] = sphere(25);
    aSphere = 1.5*[x(:)'+curState(1);y(:)'+curState(2);z(:)'+curState(3)];
    d2ActualPixels = stStereoModel.matActualProjection * [aSphere; ones(1,size(aSphere,2))];
    d2DesirePixels = stStereoModel.matDesireProjection * [aSphere; ones(1,size(aSphere,2))];
%     d3Position = curState(1:3); d3Orientation = [sin(curState(4))*cos(curState(5)); sin(curState(4))*sin(curState(5)); cos(curState(4))];
%     bodyShape = GenerateBodyShape(GenerateBodyProfile(100), d3Orientation, d3Position, 100);
%     d2ActualPixels = stStereoModel.matActualProjection * [bodyShape.points; ones(1,size(bodyShape.points,2))];
%     d2DesirePixels = stStereoModel.matDesireProjection * [bodyShape.points; ones(1,size(bodyShape.points,2))];
    d2ActualPixels(1:2, :) = d2ActualPixels(1:2, :) ./ repmat(d2ActualPixels(3, :), 2, 1); d2ActualPixels = floor(d2ActualPixels(1:2,:)); d2ActualPixels = unique(d2ActualPixels','rows')';
    d2DesirePixels(1:2, :) = d2DesirePixels(1:2, :) ./ repmat(d2DesirePixels(3, :), 2, 1); d2DesirePixels = floor(d2DesirePixels(1:2,:)); d2DesirePixels = unique(d2DesirePixels','rows')';
    
%     curActualPt = stStereoModel.matActualProjection * [curState(1:3);1]; curActualPt = curActualPt /curActualPt(3);
%     curDesirePt = stStereoModel.matDesireProjection * [curState(1:3);1]; curDesirePt = curDesirePt /curDesirePt(3);
    %-------------------------------------------------
    subplot(121); hold on; imshow(subActualImage);
	plot(scale*(d2ActualPixels(1,:)-curActualPixel(1)+subSize), scale*(d2ActualPixels(2,:)-curActualPixel(2)+subSize),'.y');
% 	plot(scale*(curActualPt(1)-curActualPixel(1)+subSize),scale*(curActualPt(2)-curActualPixel(2)+subSize),'og');
    subplot(122); hold on; imshow(subDesireImage);
%     plot(scale*(curDesirePt(1)-curDesirePixel(1)+subSize),scale*(curDesirePt(2)-curDesirePixel(2)+subSize),'og');
	plot(scale*(d2DesirePixels(1,:)-curDesirePixel(1)+subSize), scale*(d2DesirePixels(2,:)-curDesirePixel(2)+subSize),'.y');
end

%
%
%
function DrawParticle(stStereoModel, curState, inParticles, inWeights, color)

    %--------------------------
    idxes = find(inWeights>=mean(inWeights));
    tm0Positions = inParticles(1:3, idxes); tm0Positions(4,:) = 1;
    %--------------------------
    
    tm0ActualPts = stStereoModel.matActualProjection * tm0Positions;
    tm0ActualPts(1:2, :) = tm0ActualPts(1:2, :) ./ repmat(tm0ActualPts(3, :), 2, 1);
    tm0DesirePts = stStereoModel.matDesireProjection * tm0Positions;
    tm0DesirePts(1:2, :) = tm0DesirePts(1:2, :) ./ repmat(tm0DesirePts(3, :), 2, 1);
    %--------------------------------------------------------
	subplot(121); %hold on; plot(tm0ActualPts(1,:), tm0ActualPts(2,:), '.', 'color',color);
    curActualPixel = stStereoModel.matActualProjection * [curState(1:3);1]; curActualPixel = curActualPixel /curActualPixel(3);
    plot(curActualPixel(1),curActualPixel(2),'og');
    %------------------------------------------------
    subplot(122); %hold on; plot(tm0DesirePts(1,:), tm0DesirePts(2,:), '.', 'color',color);
    curDesirePixel = stStereoModel.matDesireProjection * [curState(1:3);1]; curDesirePixel = curDesirePixel /curDesirePixel(3);
    plot(curDesirePixel(1),curDesirePixel(2),'og');
end