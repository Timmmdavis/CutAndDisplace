function [ Tnn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FaceNormalVector,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz )
%CalculateNormalAndShearTractions3d Calculates normal and shear tractions
%from input direction cosines and stress tensors
%   FaceNormalVector is a 3*n vector of direction cosines: cosAx,cosAy,cosAz
%   Pxx etc are stress tensors
%   This func uses these to calc stress in normal (Tnn) strike slip (Tss) and
%   dip slip directions (Tds)

%   Copyright 2017, Tim Davis, The University of Aberdeen

    [ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector ) ;  
    %Splitting the face normal vector into its direction cosines. Note these are kept as radians not degrees. 
    CosAx=FaceNormalVector(:,1); 
    CosAy=FaceNormalVector(:,2);
    CosAz=FaceNormalVector(:,3);
    %Calculating the normal stresses on the planes
    [ Tnn ] = CalculateNormalTraction3d( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz );
    %Turning these stress vectors into traction components 
    [ Tx,Ty,Tz ] = TractionVectorCartesianComponents3d(  Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz );
    %%%%%%%%%%%%%%%%%%%%%StrikeSlipTraction calculation
    [ Tss ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,StrikeSlipCosine );
    %%%%%%%%%%%%%%%%%%%%%DipSlipStress
    [ Tds ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,DipSlipCosine );


end

