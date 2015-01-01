close all
clc;

%%
clear all
Datasets = [  5,2,1001;
             10,2,1001;
             20,2,1001;
             40,2,1001;
             80,2,1001];  
ds = 2;
DatasetPath = sprintf('simu%03d',Datasets(ds, 1));
load (sprintf('simu%03d\\stereoModel.mat', Datasets(ds, 1)));
StartImageSeq = Datasets(ds,2);
isdebug = 1;
ismultijob = 0;

%%
%
disp('tracking performance of trackers with the help of orientation');
clear trackers trajectories;
%
load(sprintf('d3TrackingData_%03d_t%04d_1.mat',Datasets(ds,1),Datasets(ds,3)-Datasets(ds,2)+1));
load(sprintf('simu%03d\\trajectories%03d.mat', Datasets(ds,1),Datasets(ds,1)));

for t = Datasets(ds,2) : Datasets(ds,3)
    grdLocations{t} = [];
    for o = 1 : Datasets(ds,1)
        grdLocations{t} = [grdLocations{t}, trajectories(o).pts(:, t-trajectories(o).start+1)];
    end
end

dis = 0; cnt = 0;
for o = 1 : length(trackers)
    obj = trackers(o);
    
    for t = obj.start : obj.end
        diss = distance(obj.states(1:3, t-obj.start+1), grdLocations{t+1});
        [m, i] = min(diss);
        if ( m<= stTrackingParameter.lenBody )
            trackers(o).grd(t-obj.start+1) = i;
            
            dis = dis + m; cnt = cnt + 1;
            trajectories(i).tracker(t+1) = o;
        else
            trackers(o).grd(t-obj.start+1) = 0;
            trajectories(i).tracker(t+1) = 0;
        end
    end
end
for i = 1 : length(trajectories)
    if ( length(trajectories(i).tracker) < Datasets(ds,3) )
        trajectories(i).tracker(length(trajectories(i).tracker)+1 : Datasets(ds,3)) = 0;
    end
end
% figure(1); hold on;
% colors = jet(length(trackers)+1);
% for o = 1 : length(trackers)
%     plot([trackers(o).start:trackers(o).end],trackers(o).grd, '-','color',colors(o,:), 'linewidth', 1);
% end

%------------------------------------------------------------------------
MOTP = dis/cnt

misses = 0;
for i = 1 : length(trajectories)
    miss = find(trajectories(i).tracker(Datasets(ds,2):Datasets(ds,3))==0);
    
    misses = misses + length(miss);
end

fpostives = 0; mismatches = 0;
for o = 1 : length(trackers)
    fp = find(trackers(o).grd==0);
    
    fpostives = fpostives + length(fp);
    
    mme = unique(trackers(o).grd); mme = setdiff(mme, 0);
    mismatches = mismatches + numel(mme) - 1;
end

MOTA = 1 - (misses+fpostives+mismatches) / (Datasets(ds,1)*(Datasets(ds,3)-Datasets(ds,2)+1))

%%
%
disp('tracking performance of trackers without the help of orientation');
clear trackers trajectories;
%
load(sprintf('d3TrackingData_%03d_t%04d_0.mat',Datasets(ds,1),Datasets(ds,3)-Datasets(ds,2)+1));
load(sprintf('simu%03d\\trajectories%03d.mat', Datasets(ds,1),Datasets(ds,1)));



for t = Datasets(ds,2) : Datasets(ds,3)
    grdLocations{t} = [];
    for o = 1 : Datasets(ds,1)
        grdLocations{t} = [grdLocations{t}, trajectories(o).pts(:, t-trajectories(o).start+1)];
    end
end

dis = 0; cnt = 0;
for o = 1 : length(trackers)
    obj = trackers(o);
    
    for t = obj.start : obj.end
        diss = distance(obj.states(1:3, t-obj.start+1), grdLocations{t+1});
        [m, i] = min(diss);
        if ( m<= stTrackingParameter.lenBody )
            trackers(o).grd(t-obj.start+1) = i;
            
            dis = dis + m; cnt = cnt + 1;
            trajectories(i).tracker(t+1) = o;
        else
            trackers(o).grd(t-obj.start+1) = 0;
            trajectories(i).tracker(t+1) = 0;
        end
    end
end
% figure(1); hold on;
% colors = jet(length(trackers)+1);
% for o = 1 : length(trackers)
%     plot([trackers(o).start:trackers(o).end],trackers(o).grd, '-','color',colors(o,:), 'linewidth', 1);
% end

%%
MOTP = dis/cnt

misses = 0;
for i = 1 : length(trajectories)
    miss = find(trajectories(i).tracker(Datasets(ds,2):Datasets(ds,3))==0);
    
    misses = misses + length(miss);
end

fpostives = 0; mismatches = 0;
for o = 1 : length(trackers)
    fp = find(trackers(o).grd==0);
    
    fpostives = fpostives + length(fp);
    
    mme = unique(trackers(o).grd); mme = setdiff(mme, 0);
    mismatches = mismatches + numel(mme) - 1;
end

MOTA = 1 - (misses+fpostives+mismatches) / (Datasets(ds,1)*(Datasets(ds,3)-Datasets(ds,2)+1))



%%
motp1 = [1.11,1.23,1.38,1.39,1.34];
mota1 = [0.99,0.98,0.96,0.94,0.94];
motp2 = [1.32,1.32,1.31,1.51,1.62];
mota2 = [0.91,0.90,0.90,0.76,0.51];

s = [5,10,20,40,80];

figure1 = figure(2); clf; set(gcf,'Position',[500,200,1400,600]);  
subplot(121); hold on; grid on;
l1 = plot(s,motp1, '-or', 'linewidth', 2);
l2 = plot(s,motp2, '-^b', 'linewidth', 2);
hLegend1 = legend([l1,l2],'3D-LOT','3D-LT'); set(hLegend1, 'fontsize', 24);
xlabel('# targets', 'fontsize', 24); ylabel('MOTP (mm)', 'fontsize', 24);
text(220,0.1,'(a)','fontsize',30);
axis([0,85,0,2]); 
%set(gca, 'Units', 'normalized', 'Position', [0.1 0.1 0.87 0.89]);  
set(gca, 'Units', 'normalized', 'Position', [0.05 0.1 0.44 0.89]);  

%figure(4); clf; set(gcf,'Position',[500,200,700,600]);  hold on; grid on;
subplot(122); hold on; grid on;
l1 = plot(s,100*mota1, '-or', 'linewidth', 2);
l2 = plot(s,100*mota2, '-^b', 'linewidth', 2);
hLegend2 = legend([l1,l2],'3D-LOT','3D-LT'); set(hLegend2, 'fontsize', 24);
xlabel('# targets', 'fontsize', 24); ylabel('MOTA (%)', 'fontsize', 24);
text(220,5,'(b)','fontsize',30);
axis([0,85,0,100]); 
%set(gca, 'Units', 'normalized', 'Position', [0.1 0.1 0.87 0.89]);  
set(gca, 'Units', 'normalized', 'Position', [0.55 0.1 0.44 0.89]);  

%%
set(hLegend1,...
    'Position',[0.33882380952381 0.319722222222225 0.1264 0.136666666666667],...
    'FontSize',24);
set(hLegend2,...
    'Position',[0.838823809523809 0.319722222222224 0.1264 0.136666666666667],...
    'FontSize',24);

% Create textbox
annotation(figure1,'textbox',...
    [0.411428571428571 0.0973333333333335 0.055 0.111666666666667],...
    'String',{'(a)'},...
    'FontSize',36,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.91 0.0973333333333339 0.055 0.111666666666667],'String',{'(b)'},...
    'FontSize',36,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'LineStyle','none');



