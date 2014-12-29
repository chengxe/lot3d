function [cross, lambda] = LineCrossEllipse(ellipse, line)
% linePoint = line.p0 + lambda*line.d
% ellipse.e, ellipse.center
    
% (linePoint-ellipse.center)' * ellipse.e * (linePoint-ellipse.center) = 1
    

    lp1 = line.p0(1); lp2 = line.p0(2); d1 = line.d(1); d2 = line.d(2); 
    a = ellipse.e(1,1); b = ellipse.e(2,2); c = ellipse.e(1,2); ep1 = ellipse.center(1); ep2 = ellipse.center(2);
%     
%     t = sym('t');
%     root = solve('(a*(lp1+t*d1-ep1)+c*(lp2+t*d2-ep2)) * (lp1+t*d1-ep1) + (c*(lp1+t*d1-ep1)+b*(lp2+t*d2-ep2)) * (lp2+t*d2-ep2) - 1 = 0', 't');
%     
%    lambda = eval(root);
    lambda = [((c^2*d1^2*ep2^2 - 2*c^2*d1^2*ep2*lp2 + c^2*d1^2*lp2^2 - 2*c^2*d1*d2*ep1*ep2 + 2*c^2*d1*d2*ep1*lp2 + 2*c^2*d1*d2*ep2*lp1 - 2*c^2*d1*d2*lp1*lp2 + c^2*d2^2*ep1^2 - 2*c^2*d2^2*ep1*lp1 + c^2*d2^2*lp1^2 + 2*c*d1*d2 - a*b*d1^2*ep2^2 + 2*a*b*d1^2*ep2*lp2 - a*b*d1^2*lp2^2 + a*d1^2 + 2*a*b*d1*d2*ep1*ep2 - 2*a*b*d1*d2*ep1*lp2 - 2*a*b*d1*d2*ep2*lp1 + 2*a*b*d1*d2*lp1*lp2 - a*b*d2^2*ep1^2 + 2*a*b*d2^2*ep1*lp1 - a*b*d2^2*lp1^2 + b*d2^2)^(1/2) + a*d1*ep1 + b*d2*ep2 + c*d1*ep2 + c*d2*ep1 - a*d1*lp1 - b*d2*lp2 - c*d1*lp2 - c*d2*lp1)/(a*d1^2 + 2*c*d1*d2 + b*d2^2);
              -((c^2*d1^2*ep2^2 - 2*c^2*d1^2*ep2*lp2 + c^2*d1^2*lp2^2 - 2*c^2*d1*d2*ep1*ep2 + 2*c^2*d1*d2*ep1*lp2 + 2*c^2*d1*d2*ep2*lp1 - 2*c^2*d1*d2*lp1*lp2 + c^2*d2^2*ep1^2 - 2*c^2*d2^2*ep1*lp1 + c^2*d2^2*lp1^2 + 2*c*d1*d2 - a*b*d1^2*ep2^2 + 2*a*b*d1^2*ep2*lp2 - a*b*d1^2*lp2^2 + a*d1^2 + 2*a*b*d1*d2*ep1*ep2 - 2*a*b*d1*d2*ep1*lp2 - 2*a*b*d1*d2*ep2*lp1 + 2*a*b*d1*d2*lp1*lp2 - a*b*d2^2*ep1^2 + 2*a*b*d2^2*ep1*lp1 - a*b*d2^2*lp1^2 + b*d2^2)^(1/2) - a*d1*ep1 - b*d2*ep2 - c*d1*ep2 - c*d2*ep1 + a*d1*lp1 + b*d2*lp2 + c*d1*lp2 + c*d2*lp1)/(a*d1^2 + 2*c*d1*d2 + b*d2^2)];

   cross = [line.p0 + lambda(1)*line.d, line.p0 + lambda(2)*line.d];
   
end

%---------------------------------------------------------
%a*x^2+b*x*y+c*y^2+d*x+e*y+1=0;
%exc=ellipse.center(1), eyc=ellipse.center(2)
%t2th = tan(2*ellipse.theta); a=ellipse.radii(2); b=ellipse.radii(1);
% sym a b c d e
% solution = solve('exc=(b*e-2*c*d) / (4*a*d-b^2)','eyc=(b*d-2*a*e) / (4*a*d-b^2)',...
%                  't2th = b / (a-c)','a^2 = 2*(a*exc^2+c*eyc^2+b*exc*eyc-1) / (a+c+((a-c)^2+b^2)^(1/2))',...
%                  'b^2 = 2*(a*exc^2+c*eyc^2+b*exc*eyc-1) / (a+c-((a-c)^2+b^2)^(1/2))','a','b','c','d','e');