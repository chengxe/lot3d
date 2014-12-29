function [newTrackers, cntTracker] = NewTargetDetection(stStereoModel, stTrackingParameter, trackers, tm1Frame, tm0Frame, t, issimu)
    
    cntTracker = 0;
    %--------------------------------------------------------------
    d3CurrentPositions = [];
    
    if ( ~isempty(trackers) )
        fliesStartTime = cat(2, trackers.start); fliesEndTime = cat(2, trackers.end);
        idx1 = find(fliesStartTime<=t); idx2 = find(fliesEndTime>=t); idxActive = intersect(idx1, idx2);     
        for n = idxActive
%                 if ( trackers(n).start + trackers(n).tm0 - 1 == t )
%                     d3CurrentPositions = [d3CurrentPositions  trackers(n).states(:,trackers(n).tm0)];
%                 end
            d3CurrentPositions = [d3CurrentPositions  trackers(n).states(:,t-trackers(n).start+1)];
        end
    end    

    for v = 1:stTrackingParameter.numCamera
        neighbors.cams{v} = InlineFindAssociatedRegion(tm1Frame.cams(v).image, tm1Frame.cams(v).regions, ...
                               tm0Frame.cams(v).image, tm0Frame.cams(v).regions, d3CurrentPositions, stStereoModel.cams(v));
        if ( isempty(neighbors.cams{v}) )
            newTrackers = []; return; end
    end
   
    for n=1:length(neighbors.cams{1})
        %disp(n);
        [targets, cntTarget] = InlineRecostruction(stStereoModel, stTrackingParameter, tm1Frame, tm0Frame, neighbors.cams{1}(n), neighbors, issimu);
        if ( cntTarget == 0 )
            continue; end
        for j=1:cntTarget
            cntTracker = cntTracker + 1; newTrackers(cntTracker) = InlineCreateTracker(stTrackingParameter, targets(j), t-1, issimu); 
        end
    end
    
    if ( cntTracker > 0 )
        newTrackers = newTrackers'; 
    else
        newTrackers = [];
    end
end

%%
function regionNeighbors = InlineFindAssociatedRegion(tm1Image, tm1Regions, tm0Image, tm0Regions, activeStates, stCam)
%从tm0开始，往后匹配
    isFound = 0;
    %-----------------------------------------------------
    [Ry, Cx] = size(tm1Image);
    tm0Regions = InlineFilterOccupiedRegion(tm0Regions, activeStates, stCam, 0);
    tm1RegionsRange = cat(2, tm1Regions.d1limit);
    
    %-----------------------------
    searchRadius = 50;
    [rowGrid colGrid] = meshgrid(-searchRadius:searchRadius,-searchRadius:searchRadius);
    indexes = find(rowGrid.^2 + colGrid.^2 <= 50^2);    
    % [rBase,cBase] = ind2sub([2*searchRadius+1,2*searchRadius+1], indexes);
    rBasis = rowGrid(indexes)'; cBasis = colGrid(indexes)';
    %-----------------------------
    
    cntNeighbor = 1; 
     for i=1:length(tm0Regions)
         regionNeighbor.tm1 = -1; regionNeighbor.tm0 = -1;
         tm0Region = tm0Regions(i); 
         if ( tm0Region.isOccupied == 1 )
             continue; end
         try
             rgnsIndex = []; similarity = [];
             occupiedPixelList = sub2ind(size(tm0Image), ...
                                 floor(rBasis+tm0Region.centroid(2)),floor(cBasis+tm0Region.centroid(1)));
            idx1 = find(tm1RegionsRange(1,:)>=min(occupiedPixelList)); idx2 = find(tm1RegionsRange(2,:)<=max(occupiedPixelList));
            idxes = intersect(idx1, idx2);
            for idx=idxes
                overlapped = tm1Regions(idx).d1pixels( ismember(tm1Regions(idx).d1pixels, occupiedPixelList) );
                if ( 0 == isempty(overlapped) )
                    rgnsIndex = [rgnsIndex, idx];
                    intensitySimilarity = CalcIntensityNcc(tm0Image(tm0Region.d1pixels), tm1Image(tm1Regions(idx).d1pixels));
                    ellipseDistance = CalcEllipseDistance(tm0Region.ellipse, tm1Regions(idx).ellipse);
                    similarity = [similarity, intensitySimilarity * (1/ellipseDistance)];                    
                end
                
                clear overlapped
            end
            if ( ~isempty(rgnsIndex) ) 
                [~, mi] = max(similarity);
                regionNeighbor.tm0 = i;
                regionNeighbor.tm1 = rgnsIndex(mi);
            end
         catch
             warning('Index exceeds matrix dimensions. tm0Region.centroid(%s)',num2str(reshape(tm0Region.centroid,1,2)));
             continue;
         end          
         
        if ( (regionNeighbor.tm1 ~= -1) && (regionNeighbor.tm0 ~= -1) )
            regionNeighbors(cntNeighbor) = regionNeighbor;
            cntNeighbor = cntNeighbor + 1; isFound = 1;
        end
     end
     
     if ( isFound == 0 )
         regionNeighbors = []; end     
end

%%
%
%
function regions = InlineFilterOccupiedRegion(regions, activeStates, stCam, isDebug)
    
%     t=linspace(0,pi,25); p=linspace(0,2*pi,25); [theta,phi]=meshgrid(t,p);
%     x=1.5*sin(theta).*sin(phi); y=1.5*sin(theta).*cos(phi); z=1.5*cos(theta);
    [x,y,z] = sphere(25);
    unitSphere = [x(:)';y(:)';z(:)'];
    
    for n=1:length(regions)
        regions(n).isOccupied = 0; end
    
    for n=1:size(activeStates,2)
        if ( isDebug )
            color =  hsv2rgb(n/size(activeStates,2), 1, 0.9);        
        end
        state = activeStates(:,n);
        aSphere = bsxfun(@plus, 1.5*unitSphere, state(1:3));
        %----------------------------------------------------------------------------
        d2Pixels = stCam.projection * [aSphere; ones(1,size(aSphere,2))];
        d2Pixels(1:2, :) = d2Pixels(1:2, :) ./ repmat(d2Pixels(3, :), 2, 1); d2Pixels = floor(d2Pixels(1:2,:)); d2Pixels = unique(d2Pixels','rows')';
        if ( isDebug )
            plot(d2Pixels(1,:),d2Pixels(2,:),'.','color', color); 
        end
        
        try
            d1PixelsList = sub2ind(stCam.resolution, d2Pixels(2,:), d2Pixels(1, :));
        catch
            warning('target (%s), pixels -> max [%d;%d], min [%d;%d]',num2str(state(1:3)'),max(d2Pixels(1,:)), max(d2Pixels(2,:)), min(d2Pixels(1,:)), min(d2Pixels(2,:)));
            continue;
        end
        regionLimit = cat(2, regions.d1limit);
        idx1 = find(regionLimit(1,:)>=min(d1PixelsList)); idx2 = find(regionLimit(2,:)<=max(d1PixelsList));
        idxes = intersect(idx1, idx2);
        
        for idx=idxes
            overlapped = regions(idx).d1pixels(ismember(regions(idx).d1pixels, d1PixelsList));
            if ( 0 == isempty(overlapped) )
                regions(idx).isOccupied = 1;
                if ( isDebug )
                    DrawRegion(regions(idx),'-g'); end
            end
            clear regionPixels; clear overlapped;
        end        
        clear d1PixelsList; 
    end
end

%%
%
%
function [targets, cntTarget] = InlineRecostruction(stStereoModel, stTrackingParameter, tm1Frame, tm0Frame, actualNeighbor, neighbors, issimu)
    
    targets = []; cntTarget = 0;
    %------------------------------------------------------------------
    tm0ActualRegion = tm0Frame.cams(1).regions(actualNeighbor.tm0); %tm0Gammas = tm0ActualRegion.ellipse.radii(2) / tm0ActualRegion.ellipse.radii(1);
    tm1ActualRegion = tm1Frame.cams(1).regions(actualNeighbor.tm1); %tm1Gammas = tm1ActualRegion.ellipse.radii(2) / tm1ActualRegion.ellipse.radii(1);

    %------------------------------------------------------------------
    for v = 2 : stTrackingParameter.numCamera
        inds = InlineFindCorrespondenceByEpipolar(stStereoModel, v-1, tm0ActualRegion.mean, tm0Frame.cams(v).regions, neighbors.cams{v}, 0, 3);
        temp = cat(2, neighbors.cams{v}.tm0); corRegionTm0Inds{v} = temp(inds);
        temp = cat(2, neighbors.cams{v}.tm1); corRegionTm1Inds{v} = temp(inds);
    end
    
    vDesire = 2;
    for j = 1 : length(corRegionTm0Inds{vDesire})
        %------------------------------------------------------------------
        tm0Regions = tm0ActualRegion; 
        tm0Gammas = tm0ActualRegion.ellipse.radii(2) / tm0ActualRegion.ellipse.radii(1);
        if ( tm0ActualRegion.ellipse2pixel > stTrackingParameter.threshold.merged )
            tm0Gammas = 1 / tm0ActualRegion.ellipse2pixel; end    
        tm1Regions = tm1ActualRegion; 
        tm1Gammas = tm1ActualRegion.ellipse.radii(2) / tm1ActualRegion.ellipse.radii(1);
        if ( tm1ActualRegion.ellipse2pixel > stTrackingParameter.threshold.merged )
            tm1Gammas = 1 / tm1ActualRegion.ellipse2pixel; end
        %------------------------------------------------------------------
        
    	tm0DesireRegion = tm0Frame.cams(vDesire).regions( corRegionTm0Inds{vDesire}(j) );
        tm0Regions = [tm0Regions, tm0DesireRegion]; 
        tm0Gammas = [tm0Gammas, tm0DesireRegion.ellipse.radii(2) / tm0DesireRegion.ellipse.radii(1)];
        if ( tm0DesireRegion.ellipse2pixel > stTrackingParameter.threshold.merged )
            tm0Gammas(end) = 1 / tm0DesireRegion.ellipse2pixel; end
        %------------------------------------------------------------------
        tm1DesireRegion = tm1Frame.cams(vDesire).regions( corRegionTm1Inds{vDesire}(j) );
        tm1Regions = [tm1Regions, tm1DesireRegion]; 
        tm1Gammas = [tm1Gammas, tm1DesireRegion.ellipse.radii(2) / tm1DesireRegion.ellipse.radii(1)];
        if ( tm1DesireRegion.ellipse2pixel > stTrackingParameter.threshold.merged )
            tm1Gammas(end) = 1 / tm1DesireRegion.ellipse2pixel; end        
        
        if ( 1 == InlineCheckFollowEpipolarConstraint(stStereoModel, vDesire-1, tm1ActualRegion.centroid, tm1DesireRegion.centroid, 5) )
            tm1Position3d = BinocularReconstruction(stStereoModel.cams(1), stStereoModel.cams(vDesire), tm1ActualRegion.centroid, tm1DesireRegion.centroid);
            if ( 1 == IsOnWallOrOutCube(stTrackingParameter.cube, tm1Position3d, stTrackingParameter.onoroutThreshold) )
                continue; end
            %
            % 判断创建的目标跟其他cameras里的measurements关联
            %
            okFlag = 1;
            for v1 = vDesire+1 : stTrackingParameter.numCamera
                corRegionTm1Ind(v1) = InlineCheckAssociation(stTrackingParameter, stStereoModel.cams(v1), tm1Position3d, tm1Frame.cams(v1).regions, corRegionTm1Inds{v1});
                if ( 0 == corRegionTm1Ind(v1) )
                    okFlag = 0; break; end
                
                %------------------------------------------------------------------
                tm1Regions = [tm1Regions, tm1Frame.cams(v1).regions(corRegionTm1Ind(v1))];
                tm1Gammas = [tm1Gammas, tm1Frame.cams(v1).regions(corRegionTm1Ind(v1)).ellipse.radii(2) / tm1Frame.cams(v1).regions(corRegionTm1Ind(v1)).ellipse.radii(1)];
                if ( tm1Frame.cams(v1).regions(corRegionTm1Ind(v1)).ellipse2pixel > stTrackingParameter.threshold.merged )
                    tm1Gammas(end) = 1 / tm1Frame.cams(v1).regions(corRegionTm1Ind(v1)).ellipse2pixel; end
                
            end
            if (okFlag == 0)
                continue; end
            %--------------------------------------------------------------
            tm0Position3d = BinocularReconstruction(stStereoModel.cams(1), stStereoModel.cams(vDesire), tm0ActualRegion.centroid, tm0DesireRegion.centroid);
            if ( 1 == IsOnWallOrOutCube(stTrackingParameter.cube, tm0Position3d, stTrackingParameter.onoroutThreshold) )
                continue; end
            %
            % 判断创建的目标跟其他cameras里的measurements关联
            %
            okFlag = 1;
            for v1 = vDesire+1 : stTrackingParameter.numCamera
                corRegionTm0Ind(v1) = InlineCheckAssociation(stTrackingParameter, stStereoModel.cams(v1), tm0Position3d, tm0Frame.cams(v1).regions, corRegionTm0Inds{v1});
                if ( 0 == corRegionTm0Ind(v1) )
                    okFlag = 0; break; end
                
                %------------------------------------------------------------------
                tm0Regions = [tm0Regions, tm0Frame.cams(v1).regions(corRegionTm0Ind(v1))];
                tm0Gammas = [tm0Gammas, tm0Frame.cams(v1).regions(corRegionTm0Ind(v1)).ellipse.radii(2) / tm0Frame.cams(v1).regions(corRegionTm0Ind(v1)).ellipse.radii(1)];
                if ( tm0Frame.cams(v1).regions(corRegionTm0Ind(v1)).ellipse2pixel > stTrackingParameter.threshold.merged )
                    tm0Gammas(end) = 1 / tm0Frame.cams(v1).regions(corRegionTm0Ind(v1)).ellipse2pixel; end                
            end
            if (okFlag == 0)
                continue; end
                        
            %
            % gammas for computing orientation
            %
            [~, ix] = sort(tm1Gammas, 'descend'); ix = ix(1:2); ix = sort(ix, 'ascend');
            tm1Orientation = ReconstructOrientation(stStereoModel.cams(ix(1)), stStereoModel.cams(ix(2)),tm1Regions(ix(1)).ellipse,tm1Regions(ix(2)).ellipse, issimu);
            [~, ix] = sort(tm0Gammas, 'descend'); ix = ix(1:2); ix = sort(ix, 'ascend');
            tm0Orientation = ReconstructOrientation(stStereoModel.cams(ix(1)), stStereoModel.cams(ix(2)),tm0Regions(ix(1)).ellipse,tm0Regions(ix(2)).ellipse, issimu);

            %
            %
            %
            cntTarget = cntTarget + 1;
            targets(cntTarget).d3Positions = [tm1Position3d, tm0Position3d]; 
            targets(cntTarget).d3Orientations = [tm1Orientation, tm0Orientation];
            for v=1:stTrackingParameter.numCamera
                targets(cntTarget).cams(v).ellipse = tm0Regions(v).ellipse;
                targets(cntTarget).cams(v).intensities = GetPixelIntensity(tm0Frame.cams(v).image, tm0Regions(v).pixels);
            end
        end
    end

end

%%
%
%
function [correspondenceIdxes] = InlineFindCorrespondenceByEpipolar(stStereoModel, fmi, actualPt, desireRegions, desireNeighbors, which, error)
%
% 点P(X,Y)到直线Ax+By+C=0的距离为 |AX+BY+C|  除以 根号下(A^2+B^2)
%
% input:
%   actualPt - [x, y], a point position in actual image.
%   stStereoModel - structure, parameters after stereo calibration
%   desireRegions - struct array, 
%                                   .mean
%   desireNeighbors - neighbor association at desire view
%   error - torelence error.
%
%   

    desireCoefs = CalcEpipolarLineCoefs(stStereoModel, fmi, actualPt);
    
    correspondenceIdxes = [];
    for i=1:length(desireNeighbors)
        switch which
            case 0,tmRegion = desireRegions(desireNeighbors(i).tm0);
            case 1,tmRegion = desireRegions(desireNeighbors(i).tm1);
        end
        distance = ( [reshape(tmRegion.mean,1,2) 1] * desireCoefs ) / sqrt( desireCoefs(1)^2 + desireCoefs(2)^2 );
        if ( abs(distance) <= error ) 
            correspondenceIdxes = [correspondenceIdxes; i];
        end
    end
end

%%
%
%
function result = InlineCheckFollowEpipolarConstraint(stStereoModel, fmi, ptActual, ptDesire, error )
    result = 0;
    
    desireCoefs = CalcEpipolarLineCoefs(stStereoModel, fmi, ptActual);
    distance = abs( ( [reshape(ptDesire,1,2) 1] * desireCoefs ) / sqrt( desireCoefs(1)^2 + desireCoefs(2)^2 ));
    if ( distance <= error )
        result = 1;   end
end

%%
%
%
function corRegionInd = InlineCheckAssociation(stTrackingParameter, stCam, d3Position, camRegions, candidateInds)
    
    corRegionInd = 0;
    
    aSphere = bsxfun(@plus, stTrackingParameter.sphere, d3Position);
    d2Pixels = stCam.projection * [aSphere; ones(1,size(aSphere,2))];
    d2Pixels(1:2, :) = d2Pixels(1:2, :) ./ repmat(d2Pixels(3, :), 2, 1); d2Pixels = floor(d2Pixels(1:2,:)); d2Pixels = unique(d2Pixels','rows')';
    
    if ( ~isempty(find(d2Pixels(1,:)>stCam.resolution(1))) || ~isempty(find(d2Pixels(2,:)>stCam.resolution(2))) || ...
         ~isempty(find(d2Pixels(1,:)<1)) || ~isempty(find(d2Pixels(2,:)<1)) )
        return;
    end
    
    occupiedIdxes = []; overlappedRatio = [];
    try
        d1OccupiedPixels = sub2ind(stCam.resolution, d2Pixels(2,:), d2Pixels(1, :));
        regionLimit = cat(2, camRegions.d1limit);
        idx1 = find(regionLimit(1,:)>=min(d1OccupiedPixels)); idx2 = find(regionLimit(2,:)<=max(d1OccupiedPixels));
        idxes = intersect(idx1, idx2);
        for i=1:length(idxes)
            overlapped = camRegions(idxes(i)).d1pixels(ismember(camRegions(idxes(i)).d1pixels, d1OccupiedPixels));
            if ( 0 == isempty(overlapped) )
                occupiedIdxes = [occupiedIdxes, idxes(i)];
                overlappedRatio(length(occupiedIdxes)) = length(overlapped) / length(camRegions(idxes(i)).d1pixels);
            end
            clear overlapped;
        end
    catch
        warning('max [%d;%d], min [%d;%d]',max(d2Pixels(1,:)), max(d2Pixels(2,:)), min(d2Pixels(1,:)), min(d2Pixels(2,:)));
        return;
    end 
    
    temp = ismember(occupiedIdxes, candidateInds);
    a = occupiedIdxes(temp); b = overlappedRatio(temp);
    if ( ~isempty(a) )
        [~, ix] = max(b); corRegionInd = a(ix);
    end
end

%%
%
%
function aTracker = InlineCreateTracker(stTrackingParameter, aTarget, start, issimu)

    aTracker.start = start; aTracker.end = start+1;
    aTracker.tm0 = 2; aTracker.tm1 = 1; aTracker.activated = 1; aTracker.missing = 0;
    aTracker.referee.cams = aTarget.cams;
    if ( ~issimu )
        [theta1, phi1, ~] = cart2sph(aTarget.d3Orientations(1,1), aTarget.d3Orientations(3,1), -aTarget.d3Orientations(2,1));
    else
        [theta1, phi1, ~] = cart2sph(aTarget.d3Orientations(1,1), aTarget.d3Orientations(2,1), aTarget.d3Orientations(3,1));
    end
    aTracker.states(:,1) = [aTarget.d3Positions(:,1); theta1; phi1; zeros(stTrackingParameter.lenState-5, 1);];
    if ( ~issimu )
        [theta0, phi0, ~] = cart2sph(aTarget.d3Orientations(1,2), aTarget.d3Orientations(3,2), -aTarget.d3Orientations(2,2));
    else
        [theta0, phi0, ~] = cart2sph(aTarget.d3Orientations(1,2), aTarget.d3Orientations(2,2), aTarget.d3Orientations(3,2));
    end
    aTracker.states(:,2) = [aTarget.d3Positions(:,2); theta0; phi0; aTarget.d3Positions(:,1);];

    %----------------------------------------
    aTracker.particles = repmat(aTracker.states(:,2),1,stTrackingParameter.numParticle);
    if ( stTrackingParameter.noise.type == 1 )
        aTracker.particles(1:3,:) = aTracker.particles(1:3,:) + stTrackingParameter.noise.sigma(1)*(0.5-rand(3,stTrackingParameter.numParticle));
    else
        aTracker.particles(1:3,:) = aTracker.particles(1:3,:) + chol(stTrackingParameter.noise.sigma(2))*randn(3,stTrackingParameter.numParticle);
    end
	aTracker.weights = repmat(1/stTrackingParameter.numParticle,1,stTrackingParameter.numParticle);
end