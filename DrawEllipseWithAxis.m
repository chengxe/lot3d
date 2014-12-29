function DrawEllipseWithAxis(ellipse, style)
    hold on;
    d2Ptsoe = InlineGeneratePointOnEllipse(ellipse, 20);
    
    plot(d2Ptsoe(1,:), d2Ptsoe(2,:), style);
    major(:,1) = ellipse.center - ellipse.axis(:,2) * ellipse.radii(2);
    major(:,2) = ellipse.center + ellipse.axis(:,2) * ellipse.radii(2);

    minor(:,1) = ellipse.center - ellipse.axis(:,1) * ellipse.radii(1);
    minor(:,2) = ellipse.center + ellipse.axis(:,1) * ellipse.radii(1);
    
    plot(major(1,:), major(2,:), style);
    plot(minor(1,:), minor(2,:), style);
end


function [d2Points] = InlineGeneratePointOnEllipse(ellipse, nbPoint)

    th = linspace(0, 2*pi, nbPoint);
    pc = [cos(th); sin(th)];
    
    pe = sqrtm(inv(ellipse.e))*pc;
    d2Points = bsxfun(@plus, ellipse.center, pe);
end