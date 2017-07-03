function [ DnTn,DnTss,DnTds,DssTn,DssTss,DssTds,DdsTn,DdsTss,DdsTds,Dn_dx,Dn_dy,Dn_dz,Dss_dx,Dss_dy,Dss_dz,Dds_dx,Dds_dy,Dds_dz,NUM,...
StrikeSlipCosine,DipSlipCosine] = CalculateInfluenceMatrices3d(MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp)
%Calculates the influence matrices
%Creating influence matrices. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
%If any fixed displacements exist we need to create the matrices too
FD=any(abs(Fdisp))>0;

X=MidPoint(:,1);
Y=MidPoint(:,2);
Z=MidPoint(:,3);

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
if  isOctave==1;
MemoryCheckerOctave(NUM); 
elseif isOctave==0;
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
Ss=1;  %Positive =  right lat movement
Ds=0;  %Positive = 	reverse movement
Ts=0;  %Positive = 	closing movement


%creating some strings for the progress bar
StringHS='1/3 CalculatingDssInfMatrixHS';
StringFS='1/3 CalculatingDssInfMatrixFS';
%Passing to internal function at the base of the file
[Dssinfmatrix,DssDisplacementXYZ]=CreateCoeffsLoop(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD,StringHS,StringFS);

%Runs the script to create strain at every faults midpoint from a Dds
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is 'X' rows. 'X' being the
%number of elements
Ss=0;  %Positive =  right lat movement
Ds=1;  %Positive = 	reverse movement
Ts=0;  %Positive = 	closing movement

%creating some strings for the progress bar
StringHS='2/3 CalculatingDdsInfMatrixHS';
StringFS='2/3 CalculatingDdsInfMatrixFS';
%Passing to internal function at the base of the file
[Ddsinfmatrix,DdsDisplacementXYZ]=CreateCoeffsLoop(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD,StringHS,StringFS);

%Runs the script to create strain at every faults midpoint from a Dn
%magnitude of 1. (Strain inf matrix).
%Appended in list so each elements influence is in 'X' rows. 'X' being the
%number of elements. Then the next element is the next X amount of rows. 
Ss=0;  %Positive =  right lat movement
Ds=0;  %Positive = 	normal movement
Ts=1;  %Positive = 	closing movement

%creating some strings for the progress bar
StringHS='3/3 CalculatingDnInfMatrixHS';
StringFS='3/3 CalculatingDnInfMatrixFS';
%Passing to internal function at the base of the file
[Dninfmatrix,DnDisplacementXYZ]=CreateCoeffsLoop(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD,StringHS,StringFS);

clear halfspace P1 P2 P3 Ds Ss Ts X Y Z first i last

%Not doing this as an internal function to avoid ramping up mem usage                                                    
%STRIKE SLIP Splitting the output strain influence matrix into tensors and converting this to stress tensor influences. 
%Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book. 
Dss_exx = Dssinfmatrix (:,1);                         Dssinfmatrix=Dssinfmatrix(:,2:6);
Dss_eyy = Dssinfmatrix (:,1);                         Dssinfmatrix=Dssinfmatrix(:,2:5);
Dss_ezz = Dssinfmatrix (:,1);                         Dssinfmatrix=Dssinfmatrix(:,2:4);
Dss_exy = Dssinfmatrix (:,1);                         Dssinfmatrix=Dssinfmatrix(:,2:3);
Dss_exz = Dssinfmatrix (:,1);                         Dssinfmatrix=Dssinfmatrix(:,2);
Dss_eyz = Dssinfmatrix (:,1);                         clear Dssinfmatrix 
[Sxxss,Syyss,Szzss,Sxyss,Sxzss,Syzss] = HookesLaw3dStrain2Stress(Dss_exx,Dss_eyy,Dss_ezz,Dss_exy,Dss_exz,Dss_eyz,lambda,mu);
clear Dss_exx Dss_eyy Dss_ezz Dss_exy Dss_exz Dss_eyz 
if FD==1
Dss_dx = DssDisplacementXYZ (:,1);                    DssDisplacementXYZ=DssDisplacementXYZ(:,2:3);
Dss_dy = DssDisplacementXYZ (:,1);                    DssDisplacementXYZ=DssDisplacementXYZ(:,2);
Dss_dz = DssDisplacementXYZ (:,1);                    clear DssDisplacementXYZ
end
%DIP SLIP Splitting the output strain influence matrix into tensors and converting this to stress tensor influences. 
%Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book.   
Dds_exx = Ddsinfmatrix (:,1);                         Ddsinfmatrix=Ddsinfmatrix(:,2:6);
Dds_eyy = Ddsinfmatrix (:,1);                         Ddsinfmatrix=Ddsinfmatrix(:,2:5);
Dds_ezz = Ddsinfmatrix (:,1);                         Ddsinfmatrix=Ddsinfmatrix(:,2:4);
Dds_exy = Ddsinfmatrix (:,1);                         Ddsinfmatrix=Ddsinfmatrix(:,2:3);
Dds_exz = Ddsinfmatrix (:,1);                         Ddsinfmatrix=Ddsinfmatrix(:,2);
Dds_eyz = Ddsinfmatrix (:,1);                         clear Ddsinfmatrix
[Sxxds,Syyds,Szzds,Sxyds,Sxzds,Syzds] = HookesLaw3dStrain2Stress(Dds_exx,Dds_eyy,Dds_ezz,Dds_exy,Dds_exz,Dds_eyz,lambda,mu);
clear Dds_exx Dds_eyy Dds_ezz Dds_exy Dds_exz Dds_eyz 
if FD==1
Dds_dx = DdsDisplacementXYZ (:,1);                    DdsDisplacementXYZ=DdsDisplacementXYZ(:,2:3);
Dds_dy = DdsDisplacementXYZ (:,1);                    DdsDisplacementXYZ=DdsDisplacementXYZ(:,2);
Dds_dz = DdsDisplacementXYZ (:,1);                    clear DdsDisplacementXYZ
end
%TENSILE SLIP Splitting the output strain influence matrix into tensors and converting this to stress tensor influences. 
%Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book. 
Dn_exx = Dninfmatrix (:,1);                         Dninfmatrix=Dninfmatrix(:,2:6);
Dn_eyy = Dninfmatrix (:,1);                         Dninfmatrix=Dninfmatrix(:,2:5);
Dn_ezz = Dninfmatrix (:,1);                         Dninfmatrix=Dninfmatrix(:,2:4);
Dn_exy = Dninfmatrix (:,1);                         Dninfmatrix=Dninfmatrix(:,2:3);
Dn_exz = Dninfmatrix (:,1);                         Dninfmatrix=Dninfmatrix(:,2);
Dn_eyz = Dninfmatrix (:,1);                         clear Dninfmatrix
[Sxxts,Syyts,Szzts,Sxyts,Sxzts,Syzts] = HookesLaw3dStrain2Stress(Dn_exx,Dn_eyy,Dn_ezz,Dn_exy,Dn_exz,Dn_eyz,lambda,mu);
clear Dn_exx Dn_eyy Dn_ezz Dn_exy Dn_exz Dn_eyz 
if FD==1
Dn_dx = DnDisplacementXYZ (:,1);                    DnDisplacementXYZ=DnDisplacementXYZ(:,2:3);
Dn_dy = DnDisplacementXYZ (:,1);                    DnDisplacementXYZ=DnDisplacementXYZ(:,2);
Dn_dz = DnDisplacementXYZ (:,1);                    clear DnDisplacementXYZ
end

%Reshape column vectors of stress into square matrices, then clear column vectors
dimx = NUM;
dimy = NUM;
%internal func down the bottom of this main one
[Sxxss,Syyss,Szzss,Sxyss,Sxzss,Syzss]=SquareReshapeTensors(Sxxss,Syyss,Szzss,Sxyss,Sxzss,Syzss,dimx,dimy);
[Sxxds,Syyds,Szzds,Sxyds,Sxzds,Syzds]=SquareReshapeTensors(Sxxds,Syyds,Szzds,Sxyds,Sxzds,Syzds,dimx,dimy);
[Sxxts,Syyts,Szzts,Sxyts,Sxzts,Syzts]=SquareReshapeTensors(Sxxts,Syyts,Szzts,Sxyts,Sxzts,Syzts,dimx,dimy);
if FD==1
[Dss_dx,Dss_dy,Dss_dz]=SquareReshapeDisp(Dss_dx,Dss_dy,Dss_dz,dimx,dimy);  
[Dds_dx,Dds_dy,Dds_dz]=SquareReshapeDisp(Dds_dx,Dds_dy,Dds_dz,dimx,dimy);  
[Dn_dx,Dn_dy,Dn_dz]=SquareReshapeDisp(Dn_dx,Dn_dy,Dn_dz,dimx,dimy);  
end

%Splitting the face normal vector into its direction cosines. Note these are kept as radians not degrees. 
CosAx=FaceNormalVector(:,1); 
CosAy=FaceNormalVector(:,2);
CosAz=FaceNormalVector(:,3);

% Now calculating influence normal tractions on the planes. 
%StrikeSlipDisplacement_TractionNormal
[ DssTn ] = CalculateNormalTraction3d( Sxxss,Syyss,Szzss,Sxyss,Sxzss,Syzss,CosAx,CosAy,CosAz );

%DipSlipDisplacement_TractionNormal
[ DdsTn ] = CalculateNormalTraction3d( Sxxds,Syyds,Szzds,Sxyds,Sxzds,Syzds,CosAx,CosAy,CosAz );

%TensileDisplacement_TractionNormal
[ DnTn ] = CalculateNormalTraction3d( Sxxts,Syyts,Szzts,Sxyts,Sxzts,Syzts,CosAx,CosAy,CosAz );

%Turning these stress inf matrix's to cartesian traction inf matrix's
[ DssT1x,DssT2y,DssT3z ] = TractionVectorCartesianComponents3d(Sxxss,Syyss,Szzss,Sxyss,Sxzss,Syzss,CosAx,CosAy,CosAz);
clear Sxxss Syyss Szzss Sxyss Sxzss Syzss
[ DdsT1x,DdsT2y,DdsT3z ] = TractionVectorCartesianComponents3d(Sxxds,Syyds,Szzds,Sxyds,Sxzds,Syzds,CosAx,CosAy,CosAz);
clear Sxxds Syyds Szzds Sxyds Sxzds Syzds
[ DnT1x,DnT2y,DnT3z ] = TractionVectorCartesianComponents3d(Sxxts,Syyts,Szzts,Sxyts,Sxzts,Syzts,CosAx,CosAy,CosAz);
clear Sxxts Syyts Szzts Sxyts Sxzts Syzts

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
Dss_dx = [];
Dss_dy = [];
Dss_dz = [];
	
Dds_dx = [];
Dds_dy = [];
Dds_dz = [];

Dn_dx = [];
Dn_dy = [];
Dn_dz = [];
end
end

function [infmatrix,DisplacementXYZ]=CreateCoeffsLoop(infmatrix,DisplacementXYZ,...
    NUM,X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,nu,halfspace,FD,StringHS,StringFS)
%Loop that calls the TDE functions of Mehdi Nikkhoo and fills large column
%coeff matrices
if halfspace==1
	progressbar(StringHS) % Create figure and set starting time
	for i=1:size(P1,1)
		first = (i-1)*NUM+1;
		last = i*NUM;
		infmatrix(first:last,:) = TDstrain_stressHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss,Ds,Ts,mu,lambda);
        if FD==1
        DisplacementXYZ(first:last,:) = TDdispHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss,Ds,Ts,nu);
        end
		progressbar(i/NUM) % Update figure 
	end
else 
	progressbar(StringFS) % Create figure and set starting time
	for i=1:size(P1,1)
		first = (i-1)*NUM+1;
		last = i*NUM;
		infmatrix(first:last,:) = TDstrain_stressFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss,Ds,Ts,mu,lambda);
        if FD==1
		DisplacementXYZ(first:last,:) = TDdispFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss,Ds,Ts,nu); 
        end
		progressbar(i/NUM) % Update figure 
	end
end

end


function [Sxx,Syy,Szz,Sxy,Sxz,Syz]=SquareReshapeTensors(Sxx,Syy,Szz,Sxy,Sxz,Syz,dimx,dimy)
Sxx = reshape(Sxx,dimx,dimy);
Syy = reshape(Syy,dimx,dimy);
Szz = reshape(Szz,dimx,dimy);
Sxy = reshape(Sxy,dimx,dimy);
Sxz = reshape(Sxz,dimx,dimy);
Syz = reshape(Syz,dimx,dimy);
end

function [Ux,Uy,Uz]=SquareReshapeDisp(Ux,Uy,Uz,dimx,dimy)
Ux = reshape(Ux,dimx,dimy);
Uy = reshape(Uy,dimx,dimy);
Uz = reshape(Uz,dimx,dimy);
end


