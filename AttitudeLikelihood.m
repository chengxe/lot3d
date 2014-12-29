function distance = AttitudeLikelihood(stCam, inImage, tm0Ellipse, refEllipse, isDebug)
    %-------------------------------------------------
    %move the center of projected to reference
    refEllipse.center = tm0Ellipse.center; 
    if ( isDebug )
         hold on; DrawEllipseWithAxis(refEllipse, '-r');  DrawEllipseWithAxis(tm0Ellipse, '-g')
    end
    %-------------------------------------------------
    try
        correspondences = InlineGenerateCorrespondence(tm0Ellipse, refEllipse, isDebug);
        d2Refs = cat(2, correspondences.ref); d2Prjs = cat(2, correspondences.prj);
        distance = sum(sqrt(sum((d2Prjs-d2Refs).^2)));
    catch
        distance = 1000;
    end
end


function correspondences = InlineGenerateCorrespondence(prjEllipse, refEllipse, isDebug)
    d2Ptsoe = InlineGeneratePointOnEllipse(prjEllipse, 31);  %points on ellipse
    d2PtsoeGrdt = gradient(d2Ptsoe);
    rotate2d = @(theta) [cos(theta), sin(theta); -sin(theta), cos(theta)];
    for n=1:size(d2Ptsoe, 2)-1
        d2Gradient = unit(d2PtsoeGrdt(:,n));
        normal = rotate2d(pi/2) * d2Gradient;

        line.d = normal; line.p0 = d2Ptsoe(:,n);
        ptCross = real(LineCrossEllipse(refEllipse, line));
            %---------------------------------
            if ( isDebug )
                p1 = d2Ptsoe(:,n) - 20*normal;
                p2 = d2Ptsoe(:,n) + 20*normal;
                plot([p1(1),p2(1)],[p1(2),p2(2)],'-b');
                plot(ptCross(1,1),ptCross(2,1),'ob');
            end
            %---------------------------------        
        correspondences(n).prj = d2Ptsoe(:,n);
        correspondences(n).ref = ptCross(:,1);
    end
end


function [d2Points] = InlineGeneratePointOnEllipse(ellipse, nbPoint)

    th = linspace(0, 2*pi, nbPoint);
    pc = [cos(th); sin(th)];
    
    pe = sqrtm(inv(ellipse.e))*pc;
    d2Points = bsxfun(@plus, ellipse.center, pe);
end
