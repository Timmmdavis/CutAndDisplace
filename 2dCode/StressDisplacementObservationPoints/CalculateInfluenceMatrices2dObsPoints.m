function [InfMats]= CalculateInfluenceMatrices2dObsPoints(halfspace,x,y,MidPoint,HalfLength,nu,E,LineNormalVector )
% CalculateInfluenceMatrices2dObsPoints: Calculates the influence matricies for the
%               Observation points for use in a linear equation. 
%
% usage #1:
% [InfMats] = CalculateInfluenceMatrices2dObsPoints(halfspace,x,y,MidPoints,HalfLength,nu,E,LineNormalVector )
%
% Arguments: (input)
%   x,y       - The XY locations of the observation points
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%   MidPoint  - The element midpoints in X and Y
%
% HalfLength  - An array of each each elements half length
%
%       nu    - The Poisson's ratio
%
%       E     - The Young's modulus
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
% Fdisp       - Flag telling the user if any elements are going to be fixed
%              (if this is the case displacement influence matricies are
%              required).
%
%
% Arguments: (output)
% InfMats           - Structure containing arrays:
%                   DsSxx,DsSyy,DsSxy...
%                   DnSxx,DnSyy,DnSxy 
%                   Influence matricies of much a displacement
%                   of one element (first part of name) effects the stress
%                   on a observation point (last part of name).
%
% Example usage:
%
% [InfMats]= CalculateInfluenceMatrices2dObsPoints(halfspace,x,y,MidPoint,HalfLength,nu,E,LineNormalVector )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Number of elements
NUM=numel(x);
NUM2 = numel(MidPoint(:,1));

%Creating empty array to be filled with influence coefficients 
InfMatrix = zeros(NUM*NUM2,3); 
Dispinfmatrix=[];%zeros(NUM*NUM2,2);

%Setting up shear disp coeff matrix
Ds = 1;
Dn = 0;
%Running loop and filling matrices with coefficients
%Simple if/else statement to create half space or non half space
%coefficients
[DsInfMatrix,~]=CreateCoeffsLoop2d(InfMatrix,Dispinfmatrix,...
    NUM,x,y,MidPoint,HalfLength,LineNormalVector,Ds,Dn,nu,E,halfspace,0);

%Setting up normal disp coeff matrix
Ds = 0;
Dn = 1;
[DnInfMatrix,~]=CreateCoeffsLoop2d(InfMatrix,Dispinfmatrix,...
    NUM,x,y,MidPoint,HalfLength,LineNormalVector,Ds,Dn,nu,E,halfspace,0);


%  Each influence array is now a huge 5*n column vectors with
%  Sxx,Syy,Szz,Ux,Uy
%  These are reshaped into 3 square matrices where each column is an
%  elements influence on every other element. 
dimx = NUM;
dimy = NUM2;

%External functions, extracting cols and reshaping
[ DsSxx,DsSyy,DsSxy ] = ExtractCols( DsInfMatrix );   
%[ DsUx,DsUy ] = ExtractCols( DsDisplacementXY ); 
[ DsSxx,DsSyy,DsSxy]=ReshapeData2d(dimx,dimy,DsSxx,DsSyy,DsSxy); 

[ DnSxx,DnSyy,DnSxy ] = ExtractCols( DnInfMatrix );   
%[ DnUx,DnUy ] = ExtractCols( DnDisplacementXY );    
[ DnSxx,DnSyy,DnSxy ]=ReshapeData2d(dimx,dimy,DnSxx,DnSyy,DnSxy); 

InfMats.DnSxx=DnSxx;

%Now for the disp influence matricies 
InfMats.DnSxx= DnSxx; clear DnSxx
InfMats.DnSyy= DnSyy; clear DnSyy
InfMats.DnSxy= DnSxy; clear DnSxy
% InfMats.DnUx = DnUx;  clear DnUx
% InfMats.DnUy = DnUy;  clear DnUy
InfMats.DsSxx= DsSxx; clear DsSxx
InfMats.DsSyy= DsSyy; clear DsSyy
InfMats.DsSxy= DsSxy; clear DsSxy
% InfMats.DsUx = DsUx;  clear DsUx
% InfMats.DsUy = DsUy;  clear DsUy

end




