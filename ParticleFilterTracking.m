function [trackers, numWorkers] = ParticleFilterTracking(stStereoModel, stTrackingParameter, trackers, tm0Frame, t,ismultijob, issimu)
    managerWorkers = 4;
    trackersPerWorker = 6;
    
    activatedIdxes = cat(2,trackers.activated);
    activatedIdxes = find(activatedIdxes>0);
    if ( isempty(activatedIdxes) )
        numWorkers = 0; return; end
    
    display(['number of activated objects : ',num2str(length(activatedIdxes))]);
    %---------------------------------------------------------------------------
    numWorkers = 1+floor( length(activatedIdxes) / trackersPerWorker);
    if ( numWorkers > managerWorkers)
        numWorkers = managerWorkers; end
    delta = 1+round(length(activatedIdxes)/numWorkers);

    for w=1:numWorkers-1
        trackerRange(w,:) = [(w-1)*delta+1,w*delta];
    end
    trackerRange(numWorkers,:) = [(numWorkers-1)*delta+1,length(activatedIdxes)];
    
	if ( ~ismultijob )
        %---------------------------
        w = 1;
        result1 = JobParticleFilterTracking(stStereoModel, stTrackingParameter, trackers(activatedIdxes(trackerRange(w,1):trackerRange(w,2))), activatedIdxes(trackerRange(w,1):trackerRange(w,2)), tm0Frame, t, 0, issimu);
        trackers(activatedIdxes(trackerRange(w,1):trackerRange(w,2))) = result1;
    else
        sched = findResource(); job = createJob(sched);
        %---------------------------
        for w=1:numWorkers
            tasks{w} = createTask(job, @JobParticleFilterTracking, 1, {stStereoModel, stTrackingParameter, trackers(activatedIdxes(trackerRange(w,1):trackerRange(w,2))), activatedIdxes(trackerRange(w,1):trackerRange(w,2)), tm0Frame,t, 0, issimu}); 
        end

        submit(job);
        waitForState(job, 'finished');
        results = getAllOutputArguments(job);
        
        if ( isempty(results) )
            get(tasks{1}, 'ErrorMessage')
            destroy(job); 
            error('empty result!');
        end
        %-------------------------------------------------------------------------
        for w=1:numWorkers
            if ( isempty(results{w}) == 1 )
                get(tasks{w}, 'ErrorMessage')
                destroy(job); 
                error('empty result! n=%d, object idex = (%s)',w,num2str(activatedIdxes(trackerRange(w,1):trackerRange(w,2))));
            end
            trackers(activatedIdxes(trackerRange(w,1):trackerRange(w,2))) = results{w};
        end
        destroy(job);
        clear job tasks
    end
end