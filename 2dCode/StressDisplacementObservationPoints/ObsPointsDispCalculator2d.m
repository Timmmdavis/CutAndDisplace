function [Ux,Uy]=ObsPointsDispCalculator2d...
    (X,Y,MidPoint,a,nu,E,halfspace,Ds,Dn,LineNormalVector)
% ObsPointsDispCalculator2d: Calculates the disp at observation
%               points by creating a influence matrix of the elements to
%               each obs point and using matrix multiplication to find the
%                displacements. Not as fast as
%               'CalculateDisplacementOnSurroundingPoints2d.m' but if looping and the
%               problem geometry is not changing this can be a good option.
%
% usage #1:
% [StressTTotal,StressTChg,StressTReg]...
%     =ObsPointsDispCalculator2d(X,Y,MidPoint,a,nu,E,halfspace...
%     ,Ds,Dn,LineNormalVector)
%
% Arguments: (input)
%      X & Y  - The observation points locations. 
%
%    MidPoint - The element midpoints in X and Y.  
%
%       a     - An array of each each elements half length
%
%       nu    - The Poisson's ratio
%
%       E     - The Young's modulus
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%     Dn,Ds   - The calculated or defined displacement of each element. 
%
% LineNormalVector - The direction cosines of each element [CosAx,CosAy].
%
% Arguments: (output)
% Disp    - From the displacement calculation 2*n 
%                   vector,[ux,uy]. 
%
% Example usage:
%
% N/A
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen



%Getting size of observation point array
n = numel(X);

%Creating the influence matricies
[InfMats]= CalculateInfluenceMatrices2dObsPoints(halfspace,X,Y,MidPoint,a,nu,E,LineNormalVector );
 
InfMats.DnUx = DnUx;  clear DnUx
InfMats.DnUy = DnUy;  clear DnUy

%Accumulate arrays 
Aux=[InfMats.DnUx,InfMats.DsUx];
Auy=[InfMats.DnUy,InfMats.DsUy];
clear InfMats

%Concatenate ready for equation 
A= [Aux;Auy];  
 
%Concatenate vector ready for equation 
B= [Dn;Ds];  

%Run linear equation system to find result of strains on the obs points
Disp=A*B;

%Now extract variables
[ Ux,Uy ] = ExtractArraysFromVector( Disp );

end
