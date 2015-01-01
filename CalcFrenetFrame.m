function [uT, uN, uB] = CalcFrenetFrame(linCenter, uIdx)
% calculate the frenet frame, 
%
% input:
%       linCenter - a center line, can be straight or curve
%       uIdx - which point on the line.
%
% output:
%        uT, uN, uB
%
%
    if ( uIdx == 1)
        p0 = linCenter(:,1) - (linCenter(:,2)-linCenter(:,1));
        uDiff1th = linCenter(:, uIdx+1) - p0;% uDiff1th = uDiff1th/2;
        uDiff2nd = linCenter(:, uIdx+1) + p0 - 2*linCenter(:, uIdx);
    else
        if ( uIdx == size(linCenter,2) )
            p1 = linCenter(:, end) + (linCenter(:, end)-linCenter(:, end-1));
            uDiff1th = p1 - linCenter(:, uIdx-1);% uDiff1th = uDiff1th/2;
            uDiff2nd = p1 + linCenter(:, uIdx-1) - 2*linCenter(:, uIdx);
        else
            uDiff1th = linCenter(:, uIdx+1) - linCenter(:, uIdx-1);% uDiff1th = uDiff1th/2;
            uDiff2nd = linCenter(:, uIdx+1) + linCenter(:, uIdx-1) - 2*linCenter(:, uIdx);
        end
    end
    
    uT = uDiff1th / norm(uDiff1th);
    if (norm(uDiff2nd) <= 1e-10)
        idx = find(uT~=0, 1, 'first'); 
%         t = uT(1);
%         sum=0; for n=1:3
%                     if (n~=idx) 
%                         uN(n) = uT(1)-t; t=uN(n);
%                         sum = sum+uT(n); end 
%         end
%         uN(idx) = sum/-uT(idx);
%         uN = uN';
        uN = [1;1;1]; uN(idx) = 0; uN(idx) = sum(uT(uN>0))/-uT(idx);
    else
        uN = cross(uDiff1th, uDiff2nd);
    end
    uN = uN / norm(uN);
   
    uB = -cross(uT, uN);
%     uB = cross(uDiff1th, uDiff2nd);
%     if ( norm(uB) ~= 0 )
%         uB = uB / norm(uB); end
%     
%     uN = cross(uT, uB);
end