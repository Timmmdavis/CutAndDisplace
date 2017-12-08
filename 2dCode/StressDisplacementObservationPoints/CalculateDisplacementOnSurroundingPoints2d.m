function [Ux,Uy] = CalculateDisplacementOnSurroundingPoints2d(X,Y,MidPoint,HalfLength,nu,E,halfspace...
    ,Ds,Dn,LineNormalVector)
% CalculateDisplacementOnSurroundingPoints2d: calculates the disp by superposition
%               This runs the calculated/defined slip for each element and
%               works out the induced disp in the surrounding points. The
%               disp for each element are summed to a the total disp
%               for each point within the loop.
%
%               Commented lines at the base of the script allow writing of
%               a table of the disp tensors. 
%
% usage #1:
% [Ux,Uy] = CalculateDisplacementOnSurroundingPoints2d(X,Y,MidPoint,HalfLength,nu,E,halfspace...
%   ,Ds,Dn,LineNormalVector)
%
% Arguments: (input)
%      x & y  - The observation points
%
%    MidPoint - The element midpoints in X and Y.
%
%  HalfLength - An array of each each elements half length.
%
%       nu    - The Poisson's ratio.
%
%       E     - The Young's modulus.
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space.
%
%     Dn,Ds   - The calculated or defined displacement of each element. 
%
% Arguments: (output)
%  Ux,Uy   - Displacement in X,Y at each observation point (XY input).
%
% Example usage:
%
% N/A
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

% Getting size of observation point array
 n = numel(X);
 
%Direction cosines
CosAx=LineNormalVector(:,1);
CosAy=LineNormalVector(:,2); 
Beta=atan2(CosAx,-CosAy);
 
%Setting up loop where the stresses from each element are appended to a master array
DispObs = zeros(n,2);
 
%Setting up progress bar
if halfspace==1
    progressbar('Calculating Displacement U HS') % Create figure and set starting time
else 
    progressbar('Calculating Displacement U FS') % Create figure and set starting time
end
 
NumPnts=size(MidPoint(:,1),1);
%Computing effect on elements. 
if halfspace==1
    for i = 1:NumPnts
        DispObs =DispObs + LDdispHS(X,Y,MidPoint(i,1),MidPoint(i,2),HalfLength(i),Beta(i),Ds(i),Dn(i),nu);    
        progressbar(i/NumPnts) % Update figure
    end
else
    for i = 1:NumPnts
        DispObs =DispObs + LDdispFS(X,Y,MidPoint(i,1),MidPoint(i,2),HalfLength(i),Beta(i),Ds(i),Dn(i),nu); 
        progressbar(i/NumPnts) % Update figure
    end
end

  %Spliiting the output from the calcuation into 3 column vectors. 
Ux = DispObs(:,1);
Uy = DispObs(:,2);

end
