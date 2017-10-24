function [ DnTn,DnTss,DnTds,DssTn,DssTss,DssTds,DdsTn,DdsTss,DdsTds,Dn_dx,Dn_dy,Dn_dz,Dss_dx,Dss_dy,Dss_dz,Dds_dx,Dds_dy,Dds_dz,NUM,...
StrikeSlipCosine,DipSlipCosine] = CalculateInfluenceMatrices3d(MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp)
%Calculates the influence matrices
%Creating influence matrices. 

%   Copyright 2017, Tim Davis, The University of Potsdam
%If any fixed displacements exist we need to create the matrices too
FD=any(abs(Fdisp))>0;

%Element midpoints Cartesian coordinates
[ X,Y,Z ] = ExtractCols( MidPoint );  

%Moving along normal vector a very small amount to get disp correct all
%the time. Tried with 'Eps', not quite enough as the cosines are below 1
%before multiplication.
X=X-(FaceNormalVector(:,1)*1e-12);
Y=Y-(FaceNormalVector(:,2)*1e-12);
Z=Z-(FaceNormalVector(:,3)*1e-12);

%Doing a memory check, will the inf matrices exceed the RAM and freeze the
%comp?
%First checking if in Octave or MATLAB, Octave has a different way of checking for free
%memory (uses Java)
%creating the size
NUM = size(X,1);
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
MemoryCheckerOctave(NUM); 
elseif isOctave==0
if FD==0
MemoryChecker(NUM,0); %not creating disp matrices
%disp('function MemoryChecker is off. Turn on. In line 31&34 function CalculateInfluencematrices3d')
else
MemoryChecker(NUM,1); %are creating, way bigger
%disp('function MemoryChecker is off. Turn on. In line 31&34 function CalculateInfluencematrices3d')
end
end

Stressinfmatrix = zeros(NUM*NUM,6); 
%Stressinfmatrix = zeros(NUM*NUM,6,'single'); disp('influence matrices in CalculateInfluencematrices3d are currently single precision')
if FD==1
Dispinfmatrix=zeros(NUM*NUM,3);
%Dispinfmatrix=zeros(NUM*NUM,3,'single');
else
Dispinfmatrix=[];
end

%Runs the script to create strain at every faults midpoint from a Dss
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is 'X' rows. 'X' being the
%number of elements
Dss=1;  
Dds=0;  
Dn=0; 

%creating some strings for the progress bar
StringHS='1/3 CalculatingDssInfMatrixHS';
StringFS='1/3 CalculatingDssInfMatrixFS';
%Passing to internal function at the base of the file
[Dssinfmatrix,DssDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD,StringHS,StringFS);

%Runs the script to create strain at every faults midpoint from a Dds
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is 'X' rows. 'X' being the
%number of elements
Dss=0;  
Dds=1;  
Dn=0;  

%creating some strings for the progress bar
StringHS='2/3 CalculatingDdsInfMatrixHS';
StringFS='2/3 CalculatingDdsInfMatrixFS';
%Passing to internal function at the base of the file
[Ddsinfmatrix,DdsDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD,StringHS,StringFS);

%Runs the script to create strain at every faults midpoint from a Dn
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is in 'X' rows. 'X' being the
%number of elements. Then the next element is the next X amount of rows. 
Dss=0; 
Dds=0; 
Dn=1; 

%creating some strings for the progress bar
StringHS='3/3 CalculatingDnInfMatrixHS';
StringFS='3/3 CalculatingDnInfMatrixFS';
%Passing to internal function at the base of the file
[Dninfmatrix,DnDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD,StringHS,StringFS);

clear halfspace P1 P2 P3 Ds Ss Ts X Y Z first i last Stressinfmatrix Dispinfmatrix

%STRIKE SLIP Splitting the output strain influence matrix into tensors and converting this to stress tensor influences. 
%Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book. 
[ DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz ] = ExtractCols( Dssinfmatrix );
[DssSxx,DssSyy,DssSzz,DssSxy,DssSxz,DssSyz] = HookesLaw3dStrain2Stress(DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz,lambda,mu);
clear DssExx DssEyy DssEzz DssExy DssExz DssEyz 
if FD==1
[ Dss_dx,Dss_dy,Dss_dz ] = ExtractCols( DssDisplacementXYZ );    
end
%DIP SLIP Splitting the output strain influence matrix into tensors and converting this to stress tensor influences. 
%Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book.  
[ DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz ] = ExtractCols( Ddsinfmatrix );
[DdsSxx,DdsSyy,DdsSzz,DdsSxy,DdsSxz,DdsSyz] = HookesLaw3dStrain2Stress(DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz,lambda,mu);
clear DdsExx DdsEyy DdsEzz DdsExy DdsExz DdsEyz
if FD==1
[ Dds_dx,Dds_dy,Dds_dz ] = ExtractCols( DdsDisplacementXYZ );        
end
%TENSILE SLIP Splitting the output strain influence matrix into tensors and converting this to stress tensor influences. 
%Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book. 
[ DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz ] = ExtractCols( Dninfmatrix );
[DnSxx,DnSyy,DnSzz,DnSxy,DnSxz,DnSyz] = HookesLaw3dStrain2Stress(DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz,lambda,mu);
clear DnExx DnEyy DnEzz DnExy DnExz DnEyz 
if FD==1
[ Dn_dx,Dn_dy,Dn_dz ] = ExtractCols( DnDisplacementXYZ );     
end

%Reshape column vectors of stress into square matrices, then clear column vectors
dimx = NUM;
dimy = NUM;
%external function, reshapes these to square arrays
[DssSxx,DssSyy,DssSzz,DssSxy,DssSxz,DssSyz]=ReshapeMultipleArrays(dimx,dimy,DssSxx,DssSyy,DssSzz,DssSxy,DssSxz,DssSyz);
[DdsSxx,DdsSyy,DdsSzz,DdsSxy,DdsSxz,DdsSyz]=ReshapeMultipleArrays(dimx,dimy,DdsSxx,DdsSyy,DdsSzz,DdsSxy,DdsSxz,DdsSyz);
[DnSxx,DnSyy,DnSzz,DnSxy,DnSxz,DnSyz]=ReshapeMultipleArrays(dimx,dimy,DnSxx,DnSyy,DnSzz,DnSxy,DnSxz,DnSyz);
if FD==1
[Dss_dx,Dss_dy,Dss_dz]=ReshapeMultipleArrays(dimx,dimy,Dss_dx,Dss_dy,Dss_dz); 
[Dds_dx,Dds_dy,Dds_dz]=ReshapeMultipleArrays(dimx,dimy,Dds_dx,Dds_dy,Dds_dz); 
[Dn_dx,Dn_dy,Dn_dz]=ReshapeMultipleArrays(dimx,dimy,Dn_dx,Dn_dy,Dn_dz); 
end

%Splitting the face normal vector into its direction cosines. Note these are kept as radians not degrees. 
[ CosAx,CosAy,CosAz ] = ExtractCols( FaceNormalVector );    

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
[Dss_dx,Dss_dy,Dss_dz,Dds_dx,Dds_dy,Dds_dz,Dn_dx,Dn_dy,Dn_dz ] = CreateBlankVars;    
end

end









