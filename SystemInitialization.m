function stTrackingParameter = SystemInitialization(numCamera, fps, tnoise, issimu, haso)

    if ( 1 == issimu )
        load scube.mat
    else
        load fcube.mat
    end
    stTrackingParameter.cube = InlineComputeCubePlanes(cube);
    stTrackingParameter.onoroutThreshold = 5;
    
    
    stTrackingParameter.fps = fps;
    stTrackingParameter.numCamera = numCamera;
    
    stTrackingParameter = InlineModelInitial(stTrackingParameter);
    stTrackingParameter.noise.type = tnoise; %1,uniform; 2,gaussian
    if ( stTrackingParameter.noise.type == 1 )
        stTrackingParameter.numParticle = 100;
    else
        stTrackingParameter.numParticle = 500;
    end

    stTrackingParameter.threshold.retry = 5;
    stTrackingParameter.threshold.resample = 0.5;
    stTrackingParameter.threshold.missing.particle = stTrackingParameter.numParticle-2;
    stTrackingParameter.threshold.missing.target = 1;
    stTrackingParameter.threshold.merged = 1.3;
    
    [x,y,z] = sphere(25);
    stTrackingParameter.sphere = 1.5*[x(:)';y(:)';z(:)'];
    stTrackingParameter.shape = GenerateBodyShape(GenerateBodyProfile(100), [0;0;1], [0;0;0], 100);
    
    stTrackingParameter.haso = haso;
    stTrackingParameter.lenBody = 2.73;
end

%
%
%
function stTrackingParameter = InlineModelInitial(stTrackingParameter)
    %------------------------------------
    stTrackingParameter.lenState = 8;
    % position X = (x_tm0,y_tm0,z_tm0), O = (phi_tm0,theta_tm0)
    % state = (X_tm0, O_tm0, X_tm1);
    stTrackingParameter.matTransition = [ 
                                        2 0 0 0 0 -1  0  0;
                                        0 2 0 0 0  0 -1  0;
                                        0 0 2 0 0  0  0 -1;
                                        0 0 0 1 0  0  0  0;
                                        0 0 0 0 1  0  0  0;
                                        1 0 0 0 0  0  0  0;
                                        0 1 0 0 0  0  0  0;
                                        0 0 1 0 0  0  0  0]; %FLE transition
                                    
    stTrackingParameter.noise.sigma = [3,3^2];   %FLE noise
end

%
%
%
function cube = InlineComputeCubePlanes(cube)

    vertex4Planes = [1,2,3,4;4,8,5,1;1,2,6,5;5,6,7,8;8,4,3,7;7,6,2,3];
    plane1(:,:,1) = [1,4;1,2]; plane1(:,:,2) = [2,1;2,3]; plane1(:,:,3) = [3,2;3,4]; plane1(:,:,4) = [4,3;4,1];
    plane2(:,:,1) = [4,1;4,8]; plane2(:,:,2) = [8,4;8,5]; plane2(:,:,3) = [5,8;5,1]; plane2(:,:,4) = [1,5;1,4];
    plane3(:,:,1) = [1,2;1,5]; plane3(:,:,2) = [2,6;2,1]; plane3(:,:,3) = [6,5;6,2]; plane3(:,:,4) = [5,1;5,6];
    plane4(:,:,1) = [5,6;5,8]; plane4(:,:,2) = [6,7;6,5]; plane4(:,:,3) = [7,8;7,6]; plane4(:,:,4) = [8,5;8,7];
    plane5(:,:,1) = [8,7;8,4]; plane5(:,:,2) = [4,8;4,3]; plane5(:,:,3) = [3,4;3,7]; plane5(:,:,4) = [7,3;7,8];
    plane6(:,:,1) = [7,6;7,3]; plane6(:,:,2) = [6,2;6,7]; plane6(:,:,3) = [2,3;2,6]; plane6(:,:,4) = [3,7;3,2];
    the6Planes{1} = plane1; the6Planes{2} = plane2; the6Planes{3} = plane3; the6Planes{4} = plane4; the6Planes{5} = plane5; the6Planes{6} = plane6;

    for p=1:6
        plane = the6Planes{p};

        norms = [];
        for n=1:4
            vectors = plane(:,:,n);
            theNorm = cross(cube.d3Vertices(:,vectors(1,2))-cube.d3Vertices(:,vectors(1,1)), ...
                            cube.d3Vertices(:,vectors(2,2))-cube.d3Vertices(:,vectors(2,1)));  theNorm = unit(theNorm);
            norms(:,n) = theNorm;  
        end
        norms(:,5) = mean(norms,2);

        cube.d3Norms(:,p) = norms(:,5);
        cube.d3PlaneCenters(:,p) = mean(cube.d3Vertices(:, vertex4Planes(p, :)), 2);
    end

end