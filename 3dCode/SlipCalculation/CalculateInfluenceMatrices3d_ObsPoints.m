function [DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz,DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz,...
 DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz]=CalculateInfluenceMatrices3d_ObsPoints(mu,lambda...
 ,X,Y,Z,P1,P2,P3,halfspace,nu)
%CalculateInfluenceMatrices3d_ObsPoints builds influence matricies for how
%much movement of a triangle strains surronding observation points. 

%   Copyright 2017, Tim Davis, The University of Potsdam

%Would be good to have a memory check here! 


Stressinfmatrix = zeros(NUM*NUM2,6); 
%Stressinfmatrix = zeros(NUM*NUM,6,'single'); disp('influence matrices in CalculateInfluencematrices3d are currently single precision')

%Dud variable used in another func
FD=0;

%Runs the script to create strain at every faults midpoint from a Dss
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is 'X' rows. 'X' being the
%number of elements
Ss=1;  
Ds=0;  
Ts=0; 

%creating some strings for the progress bar
StringHS='1/3 CalculatingDssInfMatrixHS';
StringFS='1/3 CalculatingDssInfMatrixFS';
%Passing to internal function at the base of the file
[Dssinfmatrix]=CreateCoeffsLoop3d(Stressinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD,StringHS,StringFS);

%Runs the script to create strain at every faults midpoint from a Dds
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is 'X' rows. 'X' being the
%number of elements
Ss=0;  
Ds=1;  
Ts=0;  

%creating some strings for the progress bar
StringHS='2/3 CalculatingDdsInfMatrixHS';
StringFS='2/3 CalculatingDdsInfMatrixFS';
%Passing to internal function at the base of the file
[Ddsinfmatrix]=CreateCoeffsLoop3d(Stressinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD,StringHS,StringFS);

%Runs the script to create strain at every faults midpoint from a Dn
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is in 'X' rows. 'X' being the
%number of elements. Then the next element is the next X amount of rows. 
Ss=0; 
Ds=0; 
Ts=1; 

%creating some strings for the progress bar
StringHS='3/3 CalculatingDnInfMatrixHS';
StringFS='3/3 CalculatingDnInfMatrixFS';
%Passing to internal function at the base of the file
[Dninfmatrix]=CreateCoeffsLoop3d(Stressinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD,StringHS,StringFS);

clear halfspace P1 P2 P3 Ds Ss Ts X Y Z first i last Stressinfmatrix Dispinfmatrix

%Not doing this as an internal function to avoid ramping up mem usage                                                    
%STRIKE SLIP Splitting the output strain influence matrix into tensors
[ DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz ] = ExtractCols( Dssinfmatrix );

%DIP SLIP Splitting the output strain influence matrix into tensors  
[ DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz ] = ExtractCols( Ddsinfmatrix );

%TENSILE SLIP Splitting the output strain influence matrix into tensors  
[ DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz ] = ExtractCols( Dninfmatrix );

%Reshape column vectors of stress into square matrices, then clear column vectors
dimx = NUM;
dimy = NUM2;

%External function that reshapes these to the correct size. 
[DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz]=ReshapeMultipleArrays(dimx,dimy,DssExx,DssEyy,DssEzz,DssExy,DssExz,DssEyz);
[DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz]=ReshapeMultipleArrays(dimx,dimy,DdsExx,DdsEyy,DdsEzz,DdsExy,DdsExz,DdsEyz);
[DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz]=ReshapeMultipleArrays(dimx,dimy,DnExx,DnEyy,DnEzz,DnExy,DnExz,DnEyz);

end



