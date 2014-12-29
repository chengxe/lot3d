close all
clc;

%%
clear all
Datasets = [  5,2,1001;
             10,2,1001;
             20,2,1001;
             40,2,1001;
             80,2,1001;
            160,2,1001];  
ds = 2;
DatasetPath = sprintf('simu%03d',Datasets(ds, 1));
load (sprintf('simu%03d\\stereoModel.mat', Datasets(ds, 1)));
StartImageSeq = Datasets(ds,2);
isdebug = 1;
ismultijob = 0;

%%
% initialization
stTrackingParameter = SystemInitialization(3,100,1,1,1);

%------------------------------------------------------------------------------
disp('Initializing trackers'); tic
n = StartImageSeq;
    for v=1:stTrackingParameter.numCamera
        I1 = imread(sprintf('%s\\cam%d\\im%d%03d.jpg',DatasetPath,v,v,n)); %I1 = im2double(I1); 
        tm1Frame.cams(v).image = I1;
        tm1Frame.cams(v).regions = GetMeasurement(tm1Frame.cams(v).image);
    end
%--------------------------    
n = StartImageSeq+1;
    for v=1:stTrackingParameter.numCamera
        I1 = imread(sprintf('%s\\cam%d\\im%d%03d.jpg',DatasetPath,v,v,n)); %I1 = im2double(I1); 
        tm0Frame.cams(v).image = I1;
        tm0Frame.cams(v).regions = GetMeasurement(tm0Frame.cams(v).image);
    end    
%-------------------------
rng('default');
trackers = NewTargetDetection(stStereoModel, stTrackingParameter, [], tm1Frame, tm0Frame, 2, 1);
%     %------------------------------------------------------------------------------
%     figure(9);
%     tm0Locations = []; tm1Locations = [];
%     for i = 1 : length(trackers)
%         tm0Locations = [tm0Locations, trackers(i).states(1:3,2)];
%         tm1Locations = [tm1Locations, trackers(i).states(1:3,1)];
%     end
%     for v = 1 : stTrackingParameter.numCamera
%         subplot(2,stTrackingParameter.numCamera,v); imshow(tm0Frame.cams(v).image); hold on; 
%         for m = 1 : length(tm0Frame.cams(v).regions)
%             DrawEllipseWithAxis(tm0Frame.cams(v).regions(m).ellipse, '-g');
%         end    
%         title(['t-0, v:',num2str(v)]);
%         d2Locations = stStereoModel.cams(v).projection * [tm0Locations; ones(1, size(tm0Locations, 2))];
%         d2Locations(1:2, :) = d2Locations(1:2, :) ./ repmat(d2Locations(3, :), 2, 1); 
%         plot(d2Locations(1, :), d2Locations(2, :), 'or', 'markersize', 5);
%             
%         subplot(2,stTrackingParameter.numCamera,stTrackingParameter.numCamera+v); imshow(tm1Frame.cams(v).image); hold on; 
%         for m = 1 : length(tm1Frame.cams(v).regions)
%             DrawEllipseWithAxis(tm1Frame.cams(v).regions(m).ellipse, '-g');
%         end          
%         title(['t-1, v:',num2str(v)]);
%         d2Locations = stStereoModel.cams(v).projection * [tm1Locations; ones(1, size(tm1Locations, 2))];
%         d2Locations(1:2, :) = d2Locations(1:2, :) ./ repmat(d2Locations(3, :), 2, 1); 
%         plot(d2Locations(1, :), d2Locations(2, :), 'or', 'markersize', 5);        
%     end
%     return;
%     %------------------------------------------------------------------------------
tm1Frame = tm0Frame; clear tm0Frame
save('d3TrackingData.mat','stTrackingParameter','trackers','tm1Frame');
%-------------------------
toc    
%% 
% tracking
%
for n=StartImageSeq+2 :  Datasets(ds,3)
    t = n-StartImageSeq+1; 
    display(['tracking image ',num2str(n),' at time ',num2str(t)]); tic
    %------------------------------------------------------
    for v = 1 : stTrackingParameter.numCamera
        I1 = imread(sprintf('%s\\cam%d\\im%d%03d.jpg',DatasetPath,v,v,n));
        tm0Frame.cams(v).image = I1;
        tm0Frame.cams(v).regions = GetMeasurement(tm0Frame.cams(v).image);
        
    end   
    %------------------------------ 
    [trackers, numWorkers] = ParticleFilterTracking(stStereoModel, stTrackingParameter, trackers, tm0Frame, t,ismultijob, 1);
    stTrackingParameter.cntWorker(t) = numWorkers;
    stTrackingParameter.cntTimestep = t;
    timeCost = toc
    stTrackingParameter.time4frames(t) = timeCost;

%     %-------------------------------------------------------------
%     if ( isdebug )
%         starts = cat(2, trackers.start); ends = cat(2, trackers.end);
%         current = intersect( find(starts<=t), find(ends>=t) ); d3Locations = [];
%         for idx = current
%             d3Locations = [d3Locations, trackers(idx).states(1:3, t-trackers(idx).start + 1)];
%         end
%         
%         figure(1); 
%         for v = 1 : stTrackingParameter.numCamera
%             subplot(1,stTrackingParameter.numCamera,v); imshow(tm0Frame.cams(v).image); hold on; 
%             
%             d2Locations = stStereoModel.cams(v).projection * [d3Locations; ones(1, size(d3Locations, 2))];
%             d2Locations(1:2, :) = d2Locations(1:2, :) ./ repmat(d2Locations(3, :), 2, 1); 
%             plot(d2Locations(1, :), d2Locations(2, :), 'or', 'markersize', 5);
%             title(sprintf('v=%d, t=%03d, cost: %03.2f second', v, t, timeCost));
%         end
%         pause(0.1);
%     end
%     %-------------------------------------------------------------
    display(['detecting new target, original tracker number: ', num2str(length(trackers))]);
    newTrackers = NewTargetDetection(stStereoModel, stTrackingParameter, trackers, tm1Frame, tm0Frame, t, 1);
    if ( ~isempty(newTrackers) )
        trackers = [trackers; newTrackers];
    end
    
    tm1Frame = tm0Frame; clear tm0Frame
    %-------------------------------------------------------------
    if ( ~mod(t,100) )
        save(sprintf('d3TrackingData_%03d_t%04d.mat',Datasets(ds,1),t), 'stTrackingParameter', 'trackers'); end    
    
end