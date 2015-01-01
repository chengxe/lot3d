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
for ds = 2 : 2
    disp('orientation performance');
    clear trackers trajectories;
    
    load(sprintf('d3TrackingData_%03d_t%04d_1.mat',Datasets(ds,1),Datasets(ds,3)-Datasets(ds,2)+1));
    load(sprintf('simu%03d\\trajectories%03d.mat', Datasets(ds,1),Datasets(ds,1)));
    %
    for t = Datasets(ds,2) : Datasets(ds,3)
        grdLocations{t} = [];
        for o = 1 : Datasets(ds,1)
            grdLocations{t} = [grdLocations{t}, trajectories(o).pts(:, t-trajectories(o).start+1)];
        end
    end

    angles{ds} = [];
    for o = 1 : length(trackers)
        obj = trackers(o);

        for t = obj.start : obj.end
            diss = distance(obj.states(1:3, t-obj.start+1), grdLocations{t+1});
            [m, i] = min(diss);
            if ( m<= stTrackingParameter.lenBody )
                [d3Ori(1), d3Ori(2), d3Ori(3)] = sph2cart(obj.states(4, t-obj.start+1), obj.states(5, t-obj.start+1), 1);
                d3Grd = trajectories(i).ori(:, t+1);

                angle = acos( unit(d3Ori) * unit(d3Grd) ); angle = angle * 180 / pi;
                angles{ds} = [ angles{ds}, angle];
            end
        end
    end
end

%%
figure1 = figure(1);  clf; set(gcf,'Position',[300,10,700,600]); grid on; hold on;
bins = linspace(0, 180, 181);
hs = []; lg = []; matColor = jet ( 10+2 );

for ds = 2 : 2
    %----------------------------------------------------------------------------
    [n ,~] = hist(angles{ds}, bins);
    normalized_n = n / sum(n);
    x = bins; y = normalized_n; 
    % h1 = plot(x, y, '-r', 'linewidth',2);
    cdf = min(cumsum(normalized_n), 1); 
    h1 = plot(x, cdf, '-', 'color', matColor(1+ds*2,:), 'linewidth',2);
    hs = [hs, h1]; lg = [ lg; [Datasets(ds,1)]];

    ix = find(cdf>0.98, 1); x(ix)
    plot(x(ix), 0, '*', 'color', matColor(1+ds*2,:),  'linewidth', 2, 'markersize', 10);
end
%%
hl = legend(hs, num2str(lg)); set(hl, 'fontsize', 20);
set(gca, 'XMinorGrid', 'on'); %set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorGrid', 'on');
axis([0,180,0,1]);


set(gca, 'fontsize', 16);
xlabel('angle (degree)', 'fontsize', 20);
ylabel('probability', 'fontsize',20);
set(gca, 'Units', 'normalized', 'Position', [0.1 0.1 0.8700 0.870]);    

%%
set(hl,...
    'Position',[0.611041666666666 0.653095238095238 0.28125 0.155238095238095],...
    'FontSize',20);

% Create textbox
annotation(figure1,'textbox',...
    [0.76875 0.141857142857143 0.09625 0.0957142857142857],'String',{'(c)'},...
    'FontSize',36,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'LineStyle','none');