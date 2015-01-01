function bodyShape = GenerateBodyShape(bodyProfile, d3Orientation, d3Center, nbStep)
% generate a fly's body shape
%
%   input:
%       bodyProfile - a 2-by-nbStep profile matrix
%       d3Orientation - orientations of body center line
%       d3Center - position of a generated shape center
%       nbStep - 
%
%   output:
%       bodyShape - a shape structure, with components
%                   .lenBody - length of body
%                   .center.d - orientations of center line
%                   .center.p0 - abdoman end point of a fly's body shape
%                   .center.pc - center position of center line
%                   .center.points - 3-by-nbStep, points of center line
%                   .center.phi - phi angle of center line to z axis
%                   .center.theta - theta angle of center line to x axis
%                   .points - 3-by-n, points representation of body shape
%                   .profile - body profile
%                   .nbStep
%
    bodyShape.lenBody = 2.7295;
    
    bodyShape.center.d = d3Orientation; bodyShape.center.p0 = d3Center - (bodyShape.lenBody/2) * d3Orientation; 
    bodyShape.center.pc = d3Center;
    
    %--------------------------------------------------
    % points(:,1) ÊÇÎ²£¬ points(:,end)ÊÇÍ·
    t = linspace(0, bodyShape.lenBody, nbStep);
    bodyShape.center.points = bsxfun(@plus, bodyShape.center.p0, bodyShape.center.d*t);
    
    bodyShape.center.phi = atan(sqrt(sum(d3Orientation(1:2).^2))/d3Orientation(3)); 
    if ((d3Orientation(2) == 0) && (d3Orientation(1) == 0))
        bodyShape.center.theta = 0;
    else
        bodyShape.center.theta = atan(d3Orientation(2)/d3Orientation(1));
    end
    if ( d3Orientation(1) < 0 )
        bodyShape.center.theta = pi + bodyShape.center.theta; end
    if (bodyShape.center.phi < 0) 
        bodyShape.center.phi = pi + bodyShape.center.phi; end
    
    bodyShape.nbStep = nbStep; 
    
    bodyShape.points = [];
    
    vSteps = [8,16,32,48];
    vStepNo = 4;
    v = linspace(0,2*pi,vSteps(vStepNo));
    
    LateralRatio = 1.1; 
    DorsalRatio = 1.1;    
    u = linspace(0, bodyShape.lenBody, bodyShape.nbStep);
    [uT, uN, uB] = CalcFrenetFrame(bodyShape.center.points, 10); % since center line is not a curve.
    for uIdx=1:length(u)
%         [uT, uN, uB] = CalcFrenetFrame(bodyShape.center.points, uIdx);
        rBody = interp1(bodyProfile(1,:), bodyProfile(2,:), u(uIdx));
        uvH = bsxfun(@plus, bodyShape.center.points(:, uIdx), rBody*(LateralRatio * uN*cos(v) + DorsalRatio * uB*sin(v)));
        bodyShape.points = [bodyShape.points uvH];
    end
    
    bodyShape.profile = bodyProfile;
end