function [StressInf,DispInf]= CalculateInfluenceMatrices3d(MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp)
% CalculateInfluenceMatrices3d: Calculates the influence matricies for the
%               BEM calulation. First this calculates the amount of stress
%               on each elements midpoint from a unit slip on each element
%               (Dn,Dss,Dds). This ends up with 3 large arrays where each
%               of the columns is the effect of one element on every
%               midpoint of the fault. When reshaped each midpoint becomes
%               a seperate row in the array. These are then converted to
%               traction influence matrices using Cauchys formula, this
%               uses the normal orientation of the fault at each midpoint.
%
% usage #1:
% [ StressInf,DispInf ] = CalculateInfluenceMatrices3d(MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp)
%
% Arguments: (input)
%    MidPoint - The element midpoints (triangle). (col vec, XYZ)
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
%       mu    - Shear modulus.
%
%     lambda  - Lame's constant.
%
% FaceNormalVector - The direction cosines, CosAx (Nx), CosAy and CosAz in 
%                   a list for each element. 
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
%       nu    - The Poisson's ratio.
%
%     Fdisp   - Flag telling the user if any elements are going to be fixed
%              (if this is the case displacement influence matricies are
%              required).
%
% Arguments: (output)
% StressInf          -Structure containing:
%                     DnTn,DnTss,DnTds
%                     DssTn,DssTss,DssTds
%                     DdsTn,DdsTss,DdsTds
%                     Square influence matricies of much a displacement
%                     of one element (first part of name) effects the traction
%                     on another element (last part of name).
%
% DispInf            -Structure containing:
%                     Dn_dx,Dn_dy,Dn_dz
%                     Dss_dx,Dss_dy,Dss_dz
%                     Dds_dx,Dds_dy,Dds_dz 
%                     Square influence matricies of much a displacement
%                     of one element (first part of name) effects the
%                     displacement at the midpoint of another
%                     element (not the elements displacement itself). 
%
% Example usage:
%
% [ StressInf,DispInf] = CalculateInfluenceMatrices3d(MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp)
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%If any fixed displacements exist we need to create the matrices too
FD=any(abs(Fdisp))>0;

%Element midpoints Cartesian coordinates
[ X,Y,Z ] = ExtractCols( MidPoint );

%Var NUM
NUM=size(X,1);

%Doing a memory check, will the inf matrices exceed the RAM and freeze the
%comp?
% First checking if in Octave or MATLAB, Octave has a different way of checking for free
% memory (uses Java)
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
    MemoryCheckerOctave(NUM,6);
elseif isOctave==0
    if FD==0
        MemoryChecker(NUM,6); %not creating disp matrices
        %disp('function MemoryChecker is off.
    else
        MemoryChecker(NUM,9); %are creating, way bigger
        %disp('function MemoryChecker is off.
    end
end

%Moving along normal vector a very small amount to get disp correct all
%the time. Tried with 'Eps', not quite enough as the cosines are below 1
%before multiplication.
%Splitting the face normal vector into its direction cosines. Note these are kept as radians not degrees. 
[ CosAx,CosAy,CosAz ] = ExtractCols( FaceNormalVector );    
X=X-(CosAx*1e-12);
Y=Y-(CosAy*1e-12);
Z=Z-(CosAz*1e-12);

Stressinfmatrix = zeros(NUM^2,6); 
%Stressinfmatrix = zeros(NUM^2,6,'single'); disp('using single inf mats')
if FD==1
    Dispinfmatrix=zeros(NUM^2,3);
    %Dispinfmatrix=zeros(NUM^2,3,'single');
else
    Dispinfmatrix=[];
end

%Runs the script to create strain at every faults midpoint from a Dss
%magnitude of 1. (Strain inf matrix).
Dss=1;  
Dds=0;  
Dn=0; 
%Passing to internal function at the base of the file
[Dssinfmatrix,DssDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD);

%Runs the script to create strain at every faults midpoint from a Dds
%magnitude of 1. (Strain inf matrix).
Dss=0;  
Dds=1;  
Dn=0;  
%Passing to internal function at the base of the file
[Ddsinfmatrix,DdsDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD);

%Runs the script to create strain at every faults midpoint from a Dn
%magnitude of 1. (Strain inf matrix).
Dss=0; 
Dds=0; 
Dn=1; 
%Passing to internal function at the base of the file
[Dninfmatrix,DnDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD);

clear halfspace P1 P2 P3 Ds Ss Ts X Y Z first i last Stressinfmatrix Dispinfmatrix

%StrikeSlipParts: Splitting the output strain influence matrix into tensors and converting this to stress tensor influences. 
%Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book. 
[DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz] = ExtractCols(Dssinfmatrix);
[DssSxx,DssSyy,DssSzz,DssSxy,DssSxz,DssSyz] = HookesLaw3dStrain2Stress(DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz,lambda,mu);
clear DssExx DssEyy DssEzz DssExy DssExz DssEyz Dssinfmatrix
if FD==1
    [ DssUx,DssUy,DssUz ] = ExtractCols( DssDisplacementXYZ );    
end
clear DssDisplacementXYZ

%DipSlipParts: Splitting the output strain influence matrix into tensors and converting this to stress tensor influences. 
%Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book.  
[DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz] = ExtractCols(Ddsinfmatrix);
[DdsSxx,DdsSyy,DdsSzz,DdsSxy,DdsSxz,DdsSyz] = HookesLaw3dStrain2Stress(DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz,lambda,mu);
clear DdsExx DdsEyy DdsEzz DdsExy DdsExz DdsEyz Ddsinfmatrix
if FD==1
    [ DdsUx,DdsUy,DdsUz ] = ExtractCols( DdsDisplacementXYZ );        
end
clear DdsDisplacementXYZ

%OpeningParts: Splitting the output strain influence matrix into tensors and converting this to stress tensor influences. 
%Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book. 
[DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz] = ExtractCols(Dninfmatrix);
[DnSxx,DnSyy,DnSzz,DnSxy,DnSxz,DnSyz] = HookesLaw3dStrain2Stress(DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz,lambda,mu);
clear DnExx DnEyy DnEzz DnExy DnExz DnEyz Dninfmatrix
if FD==1
    [ DnUx,DnUy,DnUz ] = ExtractCols( DnDisplacementXYZ );     
end
clear DnDisplacementXYZ

%Reshape column vectors of stress into square matrices, then clear column vectors
dimx = NUM;
dimy = NUM;
%external function, reshapes these to square arrays
[DssSxx,DssSyy,DssSzz,DssSxy,DssSxz,DssSyz]=ReshapeData2d(dimx,dimy,DssSxx,DssSyy,DssSzz,DssSxy,DssSxz,DssSyz);
[DdsSxx,DdsSyy,DdsSzz,DdsSxy,DdsSxz,DdsSyz]=ReshapeData2d(dimx,dimy,DdsSxx,DdsSyy,DdsSzz,DdsSxy,DdsSxz,DdsSyz);
[DnSxx,DnSyy,DnSzz,DnSxy,DnSxz,DnSyz]=ReshapeData2d(dimx,dimy,DnSxx,DnSyy,DnSzz,DnSxy,DnSxz,DnSyz);
if FD==1
    [DssUx,DssUy,DssUz]=ReshapeData2d(dimx,dimy,DssUx,DssUy,DssUz); 
    [DdsUx,DdsUy,DdsUz]=ReshapeData2d(dimx,dimy,DdsUx,DdsUy,DdsUz); 
    [DnUx,DnUy,DnUz]=ReshapeData2d(dimx,dimy,DnUx,DnUy,DnUz); 
end

% Now calculating influence normal tractions on the planes. 
%StrikeSlipDisplacement_TractionNormal
[ DssTn ] = CalculateNormalTraction3d( DssSxx,DssSyy,DssSzz,DssSxy,DssSxz,DssSyz,CosAx,CosAy,CosAz );

%DipSlipDisplacement_TractionNormal
[ DdsTn ] = CalculateNormalTraction3d( DdsSxx,DdsSyy,DdsSzz,DdsSxy,DdsSxz,DdsSyz,CosAx,CosAy,CosAz );

%TensileDisplacement_TractionNormal
[ DnTn ] = CalculateNormalTraction3d( DnSxx,DnSyy,DnSzz,DnSxy,DnSxz,DnSyz,CosAx,CosAy,CosAz );

%Turning these stress inf matrix's to cartesian traction inf matrix's
[ DssT1x,DssT2y,DssT3z ] = TractionVectorCartesianComponents3d(DssSxx,DssSyy,DssSzz,DssSxy,DssSxz,DssSyz,CosAx,CosAy,CosAz);
clear DssSxx DssSyy DssSzz DssSxy DssSxz DssSyz
[ DdsT1x,DdsT2y,DdsT3z ] = TractionVectorCartesianComponents3d(DdsSxx,DdsSyy,DdsSzz,DdsSxy,DdsSxz,DdsSyz,CosAx,CosAy,CosAz);
clear DdsSxx DdsSyy DdsSzz DdsSxy DdsSxz DdsSyz
[ DnT1x,DnT2y,DnT3z ] = TractionVectorCartesianComponents3d(DnSxx,DnSyy,DnSzz,DnSxy,DnSxz,DnSyz,CosAx,CosAy,CosAz);
clear DnSxx DnSyy DnSzz DnSxy DnSxz DnSyz

%Calculates the directions of the dipslip and ss directions
[ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector,CosAx,CosAy,CosAz );

%StrikeSlipDisplacement_TractionStrikeSlip
[ DssTss ] = CalculateTractionInChosenDirection3d( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
%DipSlipDisplacement_TractionStrikeSlip
[ DdsTss ] = CalculateTractionInChosenDirection3d( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
%TensileDisplacement_TractionStrikeSlip
[ DnTss  ] = CalculateTractionInChosenDirection3d( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );

%StrikeSlipDisplacement_TractionDipSlip
[ DssTds ] = CalculateTractionInChosenDirection3d( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,DipSlipCosine );
%DipSlipDisplacement_TractionDipSlip
[ DdsTds ] = CalculateTractionInChosenDirection3d( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,DipSlipCosine );
%TensileDisplacement_TractionDipSlip  
[ DnTds  ] = CalculateTractionInChosenDirection3d( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,DipSlipCosine );


clear DssT1x DssT2y DssT3z DdsT1x DdsT2y DdsT3z DnT1x DnT2y DnT3z one two three
clear CosAx CosAy CosAz 

if FD==0
    [DssUx,DssUy,DssUz,DdsUx,DdsUy,DdsUz,DnUx,DnUy,DnUz ] = CreateBlankVars;    
end

%Now putting stress influence matricies inside a structure
StressInf.DnTn = DnTn;  clear DnTn
StressInf.DnTss= DnTss; clear DnTss
StressInf.DnTds= DnTds; clear DnTds
StressInf.DssTn= DssTn; clear DssTn
StressInf.DssTss=DssTss;clear DssTss
StressInf.DssTds=DssTds;clear DssTds
StressInf.DdsTn= DdsTn; clear DdsTn
StressInf.DdsTss=DdsTss;clear DdsTss
StressInf.DdsTds=DdsTds;clear DdsTds
%Now for the disp influence matricies 
DispInf.DnUx= DnUx; clear DnUx
DispInf.DnUy= DnUy; clear DnUy
DispInf.DnUz= DnUz; clear DnUz
DispInf.DssUx=DssUx;clear DssUx
DispInf.DssUy=DssUy;clear DssUy
DispInf.DssUz=DssUz;clear DssUz
DispInf.DdsUx=DdsUx;clear DdsUx
DispInf.DdsUy=DdsUy;clear DdsUy
DispInf.DdsUz=DdsUz;clear DdsUz


end









