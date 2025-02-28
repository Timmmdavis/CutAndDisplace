function [Ux,Uy,Uz]...
    =ObsPointsDispCalculator3d(Dn,Dss,Dds,mu,lambda,X,Y,Z,P1,P2,P3,halfspace,nu)
% ObsPointsStressCalculator3d: Calculates the strain on points due a fault
%                              surface with known slip. This is slower than
%                              function
%                              'CalculateStressOnSurroundingPoints' but has
%                              the advantage you can store these influence
%                              matricies if you are running a loop where
%                              both the surface and obs points geometry
%                              stay the same.
%
% usage #1:
% [Ux,Uy,Uz]...
%     =ObsPointsDispCalculator3d(Dn,Dss,Dds,mu,lambda,X,Y,Z,P1,P2,P3,halfspace,nu)
%
% Arguments: (input)
%    X,Y & Z  - The observation points
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
%  Dss,Dds,Dn - Vectors that describe how much the elements displace in the
%               normal (Dn) and strike slip (Dss) and dipslip (Dds)
%               directions on the elements.
%
%       mu    - Shear modulus.
%
%     lambda  - Lame's constant.
%
%       nu    - The Poisson's ratio.
%
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
% Arguments: (output)
%
%
% Example usage:
%
% N/A
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

n = numel(X);

%Calculate influence matricies
[InfMats]=CalculateInfluenceMatrices3dObsPoints(mu,lambda,X,Y,Z,P1,P2,P3,halfspace,nu);

%Accumulate arrays 
Aux=[InfMats.DnUx,InfMats.DssUx,InfMats.DdsUx];
Auy=[InfMats.DnUy,InfMats.DssUy,InfMats.DdsUy];
Auz=[InfMats.DnUz,InfMats.DssUz,InfMats.DdsUz];

clear InfMats


%Concatenate ready for equation 
A= [Aux;Auy;Auz];  clear Aux Auy Auz 

%Concatenate vector ready for equation 
B= [Dn;Dss;Dds];  

%Run linear equation system to find result of strains on the obs points
Disp=A*B;

%Now extract variables
[ Ux,Uy,Uz ] = ExtractArraysFromVector( Disp );

end
