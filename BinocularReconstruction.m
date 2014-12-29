function point3d = BinocularReconstruction(camActual, camDesire, pt2dActual, pt2dDesire)
%
%
%

	lineActual = Point2Matrix(pt2dActual) * camActual.projection;
    lineDesire = Point2Matrix(pt2dDesire) * camDesire.projection;
    
    A = [lineActual(:, 1:3);
         lineDesire(:, 1:3)];
    b = [-lineActual(:, 4);
         -lineDesire(:, 4)];
     
    point3d = (A\b);
end

function mat = Point2Matrix(pt2d)
    mat = [  0       -1        pt2d(2)
             1        0       -pt2d(1)
            -pt2d(2)  pt2d(1)  0];
end