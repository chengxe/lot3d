function isornot = IsOnWallOrOutCube(cube, d3Point, epsilon)
%
% cube.d3Norms(3,6)
% cube.d3PlaneCenters(3,6);
%

    isornot = 0;
    for p=1:6
        
        distance = (d3Point-cube.d3PlaneCenters(:,p))' * cube.d3Norms(:,p);
        if ( distance < epsilon )
            isornot = 1;
            return;
        end
    end
    
end