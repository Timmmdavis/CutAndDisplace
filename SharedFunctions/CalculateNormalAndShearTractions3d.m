function [ Tn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FaceNormalVector,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz )
% CalculateNormalAndShearTractions3d: Calculates normal and shear tractions
%					from input direction cosines and stress tensors. Inputs
%					can be column vectors.
%                   Positive dipslip traction faces updip, positive strike
%                   slip faces counter clockwise from the normal in XY and
%                   normal traction is positive when facing in the
%                   direction of the normal vector.
%   
% usage #1:
% [ Tn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FaceNormalVector,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz )
%
% Arguments: (input)
% FaceNormalVector  - The normal vector of the triangles of the surface, (n*3)
%				      [CosAx,CosAy,CosAz]
%
% Pxx,Pyy,Pzz 
% Pxy,Pxz,Pyz 		- The stress tensor components on this plane.
%
% Arguments: (output)
% Tn,Tds,Tss  		- 3 column vectors of the traction magnitude in the 
%				      normal (nn), dip (ds) and strike (ss) direction.
%
% Example usage 1:
%
% % Calculating tractions for a plane dipping 45 degrees facing east under
% % tension in Sxx:
% FaceNormalVector=[cosd(45),0,cosd(45)];
% Sxx=1; Syy=0; Szz=0; Sxy=0; Sxz=0; Syz=0;
% [ Tn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FaceNormalVector,Sxx,Syy,Szz,Sxy,Sxz,Syz )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Calling function to calculate directions (function)
[ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector ) ;  
%Splitting the face normal vector into its direction cosines. Note these
%are kept as radians not degrees.
CosAx=FaceNormalVector(:,1); 
CosAy=FaceNormalVector(:,2);
CosAz=FaceNormalVector(:,3);
%Calculating the normal stresses on the planes (function)
[ Tn ] = CalculateNormalTraction3d( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz );
%Turning these stress vectors into traction components (function)
[ Tx,Ty,Tz ] = TractionVectorCartesianComponents3d( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz );
%Strike slip traction calculation (function)
[ Tss ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,StrikeSlipCosine );
%Dip slip traction calculation (function)
[ Tds ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,DipSlipCosine );


end

