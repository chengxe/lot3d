close all
clc;

%%
clear all
Datasets = [94,201,4000;
            96,201,4000;
            97,201,4000; ];
dataset = 1;
StartImageSeq = Datasets(dataset,2);

DatasetPath = sprintf('e:\\datasets\\fruitfly\\20120910\\%d\\',Datasets(dataset,1));
stStereoModel = load('e:\datasets\fruitfly\20120910\calib\stStereoModel.mat');
namePrefix = 'CoreView_Camera';

isdebug = 1;
ismultijob = 0;


%%
%
% initialization
%
stTrackingParameter = SystemInitialization(3,100,2,0);

%------------------------------------------------------------------------------
disp('Initializing gaussian background model, window size=50');
clear stBackgroundModels imageSeries
NbImagesGBM = 50;
tic
%-------------------------------------------------------------------------------------
for v=1:stTrackingParameter.numCamera
    n=StartImageSeq
    I1 = imread(sprintf('%s%s',DatasetPath,sprintf('%s%d_%04d.jpg',namePrefix,v,n))); I1 = im2double(I1); 
    imageSeries = zeros([size(I1), NbImagesGBM]);
    for t=1:NbImagesGBM
        n = t+StartImageSeq-1;
        I1 = imread(sprintf('%s%s',DatasetPath,sprintf('%s%d_%04d.jpg',namePrefix,v,n))); I1 = im2double(I1); 
        imageSeries(:,:, n-StartImageSeq+1) = I1;
    end
    stBackgroundModels{v} = GBMInitial(imageSeries, Datasets(ds,3)-Datasets(ds,2)+1);
    
    clear imageSeries;
end
toc

%------------------------------------------------------------------------------
disp('Initializing trackers'); tic
n = StartImageSeq
for v=1:stTrackingParameter.numCamera    
    I1 = imread(sprintf('%s%s',DatasetPath,sprintf('%s%d_%04d.jpg',namePrefix,v,n)));
    tm1Frame.cams(v).image = I1;
    [tm1Frame.cams(v).regions, tm1Frame.cams(v).residual] = GetObjectRegion(tm1Frame.cams(v).image, stBackgroundModels{v});
    tm1Frame.cams(v).regions = ReCalcRegions(tm1Frame.cams(v).regions, tm1Frame.cams(v).image, stBackgroundModels{v});
end    
%--------------------------
n=StartImageSeq+1;
for v=1:stTrackingParameter.numCamera    
    I1 = imread(sprintf('%s%s',DatasetPath,sprintf('%s%d_%04d.jpg',namePrefix,v,n)));
    tm0Frame.cams(v).image = I1;
    [tm0Frame.cams(v).regions, tm0Frame.cams(v).residual] = GetObjectRegion(tm0Frame.cams(v).image, stBackgroundModels{v});
    tm0Frame.cams(v).regions = ReCalcRegions(tm0Frame.cams(v).regions, tm0Frame.cams(v).image, stBackgroundModels{v});
end    
%-------------------------
rng('default');
trackers = NewTargetDetection(stStereoModel, stTrackingParameter, [], tm1Frame, tm0Frame, 2);
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
for n=StartImageSeq+2 :  Datasets(dataset,3)
    t = n-StartImageSeq+1; 
    display(['tracking image ',num2str(n),' at time ',num2str(t)]); tic
    %------------------------------------------------------
    for v = 1 : stTrackingParameter.numCamera
        hSize = floor(stBackgroundModels{v}.size(3)/2);
        if ((t-hSize>=1) && (t+hSize<=stBackgroundModels{v}.size(4)) )
         clear updateImage

         un = (t-hSize+NbImagesGBM-1+StartImageSeq-1);
         updateImage = imread(sprintf('%s%s',DatasetPath,sprintf('%s%d_%04d.jpg',namePrefix,v,n)));
         stBackgroundModels{v} = GBMUpdate(im2double(updateImage), stBackgroundModels{v});
        end        
        %------------------------------------------------------
        I1 = imread(sprintf('%s%s',DatasetPath,sprintf('%s%d_%04d.jpg',namePrefix,v,n)));
        tm0Frame.cams(v).image = I1;
        [tm0Frame.cams(v).regions, tm0Frame.cams(v).residual] = GetObjectRegion(tm0Frame.cams(v).image, stBackgroundModels{v});
        tm0Frame.cams(v).regions = ReCalcRegions(tm0Frame.cams(v).regions, tm0Frame.cams(v).image, stBackgroundModels{v});
    end
    %------------------------------ 
    [trackers, numWorkers] = ParticleFilterTracking(stStereoModel, stTrackingParameter, trackers, tm0Frame, t,ismultijob, 1);
    stTrackingParameter.cntWorker(t) = numWorkers;
    stTrackingParameter.cntTimestep = t;
    timeCost = toc
    stTrackingParameter.time4frames(t) = timeCost;

    %-------------------------------------------------------------
    if ( isdebug )
        starts = cat(2, trackers.start); ends = cat(2, trackers.end);
        current = intersect( find(starts<=t), find(ends>=t) ); d3Locations = [];
        for idx = current
            d3Locations = [d3Locations, trackers(idx).states(1:3, t-trackers(idx).start + 1)];
        end
        
        figure(1); 
        for v = 1 : stTrackingParameter.numCamera
            subplot(1,stTrackingParameter.numCamera,v); imshow(tm0Frame.cams(v).image); hold on; 
            
            d2Locations = stStereoModel.cams(v).projection * [d3Locations; ones(1, size(d3Locations, 2))];
            d2Locations(1:2, :) = d2Locations(1:2, :) ./ repmat(d2Locations(3, :), 2, 1); 
            plot(d2Locations(1, :), d2Locations(2, :), 'or', 'markersize', 5);
            title(sprintf('v=%d, t=%03d, cost: %03.2f second', v, t, timeCost));
        end
        pause(0.1);
    end
    %-------------------------------------------------------------
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