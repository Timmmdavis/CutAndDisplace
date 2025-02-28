function [StrainInf]=CalculateInfluenceMatrices3dObsPoints(mu,lambda...
 ,X,Y,Z,P1,P2,P3,halfspace,nu)
% CalculateInfluenceMatrices3dObsPoints: Calculates the influence matricies for the
%               Observation points for use in a linear equation. 
%
% usage #1:
% [StrainInf]=CalculateInfluenceMatrices3dObsPoints(mu,lambda...
%  ,X,Y,Z,P1,P2,P3,halfspace,nu)
%
% Arguments: (input)
%       mu    - Shear modulus.
%
%     lambda  - Lame's constant.
%
%       nu    - The Poisson's ratio.
%
%    X,Y & Z  - The observation points
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%
% Arguments: (output)
% StrainInf           _Structure containing:
%                   D..Exx,D..Eyy,D..Ezz
%                   D..Exy,D..Exz,D..Eyz 
%                   Influence matricies of much a displacement
%                   of one element (first part of name) effects the strain
%                   on a observation point (last part of name).
%
%
% Example usage:
%
% [StrainInf]=CalculateInfluenceMatrices3dObsPoints(mu,lambda...
%  ,X,Y,Z,P1,P2,P3,halfspace,nu)
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Would be good to have a memory check here! 
%creating the size
NUM = size(X(:),1);
NUM2 = size(P1(:,1),1);

Stressinfmatrix = zeros(NUM*NUM2,6); 
Dispinfmatrix = zeros(NUM*NUM2,3,'single'); %disp('influence matrices in CalculateInfluencematrices3d are currently single precision')

%Dud variable used in another func
FD=1;

%Runs the script to create strain at every faults midpoint from a Dss
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is 'X' rows. 'X' being the
%number of elements
Ss=1;  
Ds=0;  
Ts=0; 

%Passing to internal function at the base of the file
[Dssinfmatrix,DssDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD);

%Runs the script to create strain at every faults midpoint from a Dds
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is 'X' rows. 'X' being the
%number of elements
Ss=0;  
Ds=1;  
Ts=0;  

%Passing to internal function at the base of the file
[Ddsinfmatrix,DdsDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD);

%Runs the script to create strain at every faults midpoint from a Dn
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is in 'X' rows. 'X' being the
%number of elements. Then the next element is the next X amount of rows. 
Ss=0; 
Ds=0; 
Ts=1; 

%Passing to internal function at the base of the file
[Dninfmatrix,DnDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD);

clear halfspace P1 P2 P3 Ds Ss Ts X Y Z first i last Stressinfmatrix Dispinfmatrix

%Not doing this as an internal function to avoid ramping up mem usage                                                    
%STRIKE SLIP Splitting the output strain influence matrix into tensors
[ DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz ] = ExtractCols( Dssinfmatrix );
[ DssUx,DssUy,DssUz ] = ExtractCols( DssDisplacementXYZ );
clear Dssinfmatrix DssDisplacementXYZ

%DIP SLIP Splitting the output strain influence matrix into tensors  
[ DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz ] = ExtractCols( Ddsinfmatrix );
[ DdsUx,DdsUy,DdsUz ] = ExtractCols( DdsDisplacementXYZ );
clear Ddsinfmatrix DdsDisplacementXYZ

%TENSILE SLIP Splitting the output strain influence matrix into tensors  
[ DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz ] = ExtractCols( Dninfmatrix );
[ DnUx,DnUy,DnUz ] = ExtractCols( DnDisplacementXYZ );
clear Dninfmatrix DnDisplacementXYZ

%Reshape column vectors of stress into square matrices, then clear column vectors
dimx = NUM;
dimy = NUM2;

%External function that reshapes these to the correct size. 
[DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz]=ReshapeData2d(dimx,dimy,DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz);
[DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz]=ReshapeData2d(dimx,dimy,DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz);
[DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz]=ReshapeData2d(dimx,dimy,DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz);
[DssUx,DssUy,DssUz]=ReshapeData2d(dimx,dimy,DssUx,DssUy,DssUz);
[DdsUx,DdsUy,DdsUz]=ReshapeData2d(dimx,dimy,DdsUx,DdsUy,DdsUz);
[DnUx,DnUy,DnUz]=ReshapeData2d(dimx,dimy,DnUx,DnUy,DnUz);

%Now putting stress influence matricies inside a structure
StrainInf.DnExx= DnExx; clear DnExx
StrainInf.DnEyy= DnEyy; clear DnEyy
StrainInf.DnEzz= DnEzz; clear DnEzz
StrainInf.DnExy= DnExy; clear DnExy
StrainInf.DnExz= DnExz; clear DnExz
StrainInf.DnEyz= DnEyz; clear DnEyz
StrainInf.DnUx = DnUx;  clear DnUx
StrainInf.DnUy = DnUy;  clear DnUy
StrainInf.DnUy = DnUz;  clear DnUz

StrainInf.DssExx= DssExx; clear DssExx
StrainInf.DssEyy= DssEyy; clear DssEyy
StrainInf.DssEzz= DssEzz; clear DssEzz
StrainInf.DssExy= DssExy; clear DssExy
StrainInf.DssExz= DssExz; clear DssExz
StrainInf.DssEyz= DssEyz; clear DssEyz
StrainInf.DssUx = DssUx;  clear DssUx
StrainInf.DssUy = DssUy;  clear DssUy
StrainInf.DssUy = DssUz;  clear DssUz

StrainInf.DdsExx= DdsExx; clear DdsExx
StrainInf.DdsEyy= DdsEyy; clear DdsEyy
StrainInf.DdsEzz= DdsEzz; clear DdsEzz
StrainInf.DdsExy= DdsExy; clear DdsExy
StrainInf.DdsExz= DdsExz; clear DdsExz
StrainInf.DdsEyz= DdsEyz; clear DdsEyz
StrainInf.DdsUx = DdsUx;  clear DdsUx
StrainInf.DdsUy = DdsUy;  clear DdsUy
StrainInf.DdsUy = DdsUz;  clear DdsUz

end



